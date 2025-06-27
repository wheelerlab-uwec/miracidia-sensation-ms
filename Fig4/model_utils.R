library(broom.mixed)
library(broom)
library(lme4)
library(lmerTest)

# Function to fit model with fallback strategies
fit_model <- function(feature, data, control_level = "APW") {
  data <- data %>%
    mutate(treatment = relevel(factor(treatment), ref = control_level))

  formula <- reformulate(
    "treatment + (1 | date/side/well_col/well/particle)",
    response = feature
  )

  # Strategy 1: Try original model
  tryCatch(
    {
      clean_data <- data %>%
        filter(is.finite(.data[[feature]])) %>% # Remove Inf, -Inf, NaN
        drop_na(all_of(c(
          feature,
          "treatment",
          "date",
          "side",
          "well_col",
          "well",
          "particle"
        )))

      # Check if we have enough data left
      if (nrow(clean_data) < 1000) {
        # arbitrary minimum
        stop("Insufficient data after removing infinite values")
      }

      model <- lmer(
        formula,
        data = clean_data,
        control = lmerControl(optimizer = "bobyqa")
      )

      # Check if model is singular but still usable
      if (isSingular(model)) {
        warning(paste(
          "Singular fit for",
          feature,
          "- random effects may be unreliable"
        ))
      }

      return(
        tidy(model, effects = "fixed") %>%
          filter(str_detect(term, "^treatment")) %>%
          mutate(
            feature = feature,
            comparison = paste0(term, " vs ", control_level)
          )
      )
    },
    error = function(e1) {
      # Strategy 2: Try with cleaned data
      tryCatch(
        {
          clean_data <- data %>%
            drop_na(all_of(c(
              feature,
              "treatment",
              "date",
              "side",
              "well",
              "particle"
            )))

          model <- lmer(
            formula,
            data = clean_data,
            control = lmerControl(optimizer = "bobyqa")
          )

          if (isSingular(model)) {
            warning(paste(
              "Singular fit for",
              feature,
              "after cleaning - random effects may be unreliable"
            ))
          }

          return(
            tidy(model, effects = "fixed") %>%
              filter(str_detect(term, "^treatment")) %>%
              mutate(
                feature = feature,
                comparison = paste0(term, " vs ", control_level)
              )
          )
        },
        error = function(e2) {
          # Strategy 3: Try simplified random effects structure
          tryCatch(
            {
              simple_formula <- reformulate(
                "treatment + (1 | date)",
                response = feature
              )
              clean_data <- data %>%
                drop_na(all_of(c(feature, "treatment", "date", "side")))

              model <- lmer(
                simple_formula,
                data = clean_data,
                control = lmerControl(optimizer = "bobyqa")
              )

              return(
                tidy(model, effects = "fixed") %>%
                  filter(str_detect(term, "^treatment")) %>%
                  mutate(
                    feature = feature,
                    comparison = paste0(term, " vs ", control_level),
                    note = "simplified_random_effects"
                  )
              )
            },
            error = function(e3) {
              # Strategy 4: Fall back to fixed effects only (regular lm)
              tryCatch(
                {
                  fixed_formula <- reformulate("treatment", response = feature)
                  clean_data <- data %>%
                    drop_na(all_of(c(feature, "treatment")))

                  model <- lm(fixed_formula, data = clean_data)

                  return(
                    tidy(model) %>%
                      filter(str_detect(term, "^treatment")) %>%
                      mutate(
                        feature = feature,
                        comparison = paste0(term, " vs ", control_level),
                        note = "fixed_effects_only"
                      )
                  )
                },
                error = function(e4) {
                  # Final fallback: return NA results
                  levels <- unique(data$treatment)
                  contrasts <- setdiff(levels, control_level)

                  return(tibble(
                    effect = "fixed",
                    term = paste0("treatment", contrasts),
                    estimate = NA_real_,
                    std.error = NA_real_,
                    statistic = NA_real_,
                    df = NA_real_,
                    p.value = NA_real_,
                    feature = feature,
                    comparison = paste0(
                      "treatment",
                      contrasts,
                      " vs ",
                      control_level
                    ),
                    note = "model_failed"
                  ))
                }
              )
            }
          )
        }
      )
    }
  )
}

fit_model_droplet <- function(feature, data, control_level = "APW") {
  data <- data %>%
    mutate(treatment = relevel(factor(treatment), ref = control_level))

  formula <- reformulate(
    "treatment + (1 | date/plate/well/particle)",
    response = feature
  )

  # Strategy 1: Try original model
  tryCatch(
    {
      clean_data <- data %>%
        filter(is.finite(.data[[feature]])) %>% # Remove Inf, -Inf, NaN
        drop_na(all_of(c(
          feature,
          "treatment",
          "date",
          "plate",
          "well",
          "particle"
        )))

      # Check if we have enough data left
      if (nrow(clean_data) < 1000) {
        # arbitrary minimum
        stop("Insufficient data after removing infinite values")
      }

      model <- lmer(
        formula,
        data = clean_data,
        control = lmerControl(optimizer = "bobyqa")
      )

      # Check if model is singular but still usable
      if (isSingular(model)) {
        warning(paste(
          "Singular fit for",
          feature,
          "- random effects may be unreliable"
        ))
      }

      return(
        tidy(model, effects = "fixed") %>%
          filter(str_detect(term, "^treatment")) %>%
          mutate(
            feature = feature,
            comparison = paste0(term, " vs ", control_level)
          )
      )
    },
    error = function(e1) {
      # Strategy 2: Try with cleaned data
      tryCatch(
        {
          clean_data <- data %>%
            drop_na(all_of(c(
              feature,
              "treatment",
              "date",
              "plate",
              "well",
              "particle"
            )))

          model <- lmer(
            formula,
            data = clean_data,
            control = lmerControl(optimizer = "bobyqa")
          )

          if (isSingular(model)) {
            warning(paste(
              "Singular fit for",
              feature,
              "after cleaning - random effects may be unreliable"
            ))
          }

          return(
            tidy(model, effects = "fixed") %>%
              filter(str_detect(term, "^treatment")) %>%
              mutate(
                feature = feature,
                comparison = paste0(term, " vs ", control_level)
              )
          )
        },
        error = function(e2) {
          # Strategy 3: Try simplified random effects structure
          tryCatch(
            {
              simple_formula <- reformulate(
                "treatment + (1 | date)",
                response = feature
              )
              clean_data <- data %>%
                drop_na(all_of(c(feature, "treatment", "date", "plate")))

              model <- lmer(
                simple_formula,
                data = clean_data,
                control = lmerControl(optimizer = "bobyqa")
              )

              return(
                tidy(model, effects = "fixed") %>%
                  filter(str_detect(term, "^treatment")) %>%
                  mutate(
                    feature = feature,
                    comparison = paste0(term, " vs ", control_level),
                    note = "simplified_random_effects"
                  )
              )
            },
            error = function(e3) {
              # Strategy 4: Fall back to fixed effects only (regular lm)
              tryCatch(
                {
                  fixed_formula <- reformulate("treatment", response = feature)
                  clean_data <- data %>%
                    drop_na(all_of(c(feature, "treatment")))

                  model <- lm(fixed_formula, data = clean_data)

                  return(
                    tidy(model) %>%
                      filter(str_detect(term, "^treatment")) %>%
                      mutate(
                        feature = feature,
                        comparison = paste0(term, " vs ", control_level),
                        note = "fixed_effects_only"
                      )
                  )
                },
                error = function(e4) {
                  # Final fallback: return NA results
                  levels <- unique(data$treatment)
                  contrasts <- setdiff(levels, control_level)

                  return(tibble(
                    effect = "fixed",
                    term = paste0("treatment", contrasts),
                    estimate = NA_real_,
                    std.error = NA_real_,
                    statistic = NA_real_,
                    df = NA_real_,
                    p.value = NA_real_,
                    feature = feature,
                    comparison = paste0(
                      "treatment",
                      contrasts,
                      " vs ",
                      control_level
                    ),
                    note = "model_failed"
                  ))
                }
              )
            }
          )
        }
      )
    }
  )
}
