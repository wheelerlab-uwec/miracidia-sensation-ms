library(tidyverse)
library(here)
library(tidymodels)
library(ggtext)
library(cowplot)

usethis::use_air(vscode = TRUE)

# functions ---------------------------------------------------------------

source(here("utils", "helper_functions.R"))

# import and tidy ---------------------------------------------------------

# 1.81 GB file of raw tracking/coordinate data
# NOTE: can skip this and read the summary in line 53 or uncomment and run to see raw tracks
# zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig3", "data"), files = 'skinny_tracks.rds')

# skinny <- read_rds(here(
#   "Fig3",
#   "data",
#   "skinny_tracks.rds"
# ))

# # keep any track >5 seconds (8 FPS)
# frame_filter <- skinny |>
#   group_by(date, experiment, well_row, well_col, particle) |>
#   tally() |>
#   filter(n > 40)

# skinny_filtered <- frame_filter |>
#   select(-n) |>
#   left_join(skinny)

# chunked_data <- skinny_filtered |>
#   filter(!experiment %in% c("20241003a02rvh", "20241017a01rvh")) |>
#   group_by(date, experiment, well_row, well_col, particle) |>
#   arrange(frame, .by_group = TRUE) |>
#   group_split() |>
#   map_dfr(~ split_trajectory(.x, frame_rate = 8, chunk_duration_sec = 5))

# nest_cols <- c("date", "experiment", "well_row", "well_col", "particle")
# nested <- quick_nest(chunked_data, nest_cols)

# track_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 8,
#   chunk_size = 500,
#   pixel_to_mm = 126.5
# )

# NOTE: start here
zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig3", "data"), files = 'skinny_track_summary.rds')

track_summary <- read_rds(here(
  "Fig3",
  "data",
  "skinny_track_summary.rds"
))

# slight downward slant in the background images, so true slant is upward
# the split PNGs have a downward slant which means well_row 2 is the bottom, 1 is the middle, and 0 is the top

shifted <- track_summary |>
  unnest(c(data)) |>
  mutate(
    arena = str_extract(experiment, "a[0-9]{2}")
  ) |>
  # filter out non-moving particles (false-positives)
  filter(
    speed_mean > 632,
    net_displacement > 50,
    sd_x > 5 & sd_y > 5,
    net_displacement / path_length > 0.175
  ) |>
  mutate(
    treatment = case_when(
      arena == "a01" & well_row == 2 ~ "APW1",
      arena == "a02" & well_row == 2 ~ "APW2",
      arena == "a01" & well_row == 1 ~ "Bku",
      arena == "a01" & well_row == 0 ~ "Bga",
      arena == "a02" & well_row == 1 ~ "Bsu",
      arena == "a02" & well_row == 0 ~ "Bst",
    )
  )

###################################################
################### PLOT TRACKS ###################
###################################################

(tracks1 <- shifted |>
  filter(treatment %in% c('APW1', 'APW2', 'Bga')) |>
  mutate(
    treatment = case_when(
      str_detect(treatment, "APW") ~ "Control",
      treatment == "Bga" ~ "*B. glabrata*"
    ),
    treatment = factor(
      treatment,
      levels = c(
        "Control",
        "*B. glabrata*"
      )
    )
  ) |>
  ggplot(aes(x = x, y = y, color = frame / 8 / 60)) +
  geom_path(
    aes(group = interaction(experiment, subparticle)),
    linewidth = 0.25
  ) +
  # geom_point(size = 0.5) +
  facet_wrap(facets = vars(treatment), ncol = 2, scales = 'free_y') +
  scale_color_viridis_c(
    option = 'mako',
    breaks = c(10, 30, 50),
    guide = guide_colorbar(barwidth = 0.5, barheight = 2)
  ) +
  scale_x_continuous(
    breaks = seq(min(shifted$x), max(shifted$x), length.out = 5),
    expand = c(0, 0),
    labels = c("0", "20.25", "40.5", "60.75", "81")
  ) +
  scale_y_continuous(
    breaks = seq(min(shifted$y), max(shifted$y), length.out = 3),
    expand = c(0, 0),
    labels = c("", "", "")
  ) +
  labs(color = "Time\n(minutes)", x = ("X (mm)"), y = ("Y")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    strip.text.x = element_blank(),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank()
  ) +
  NULL)

save_plot(
  here("Fig3", "plots", "subplots", "Fig3B-bottom.pdf"),
  tracks1,
  base_width = 7.5,
  base_height = 3.5,
  device = cairo_pdf
)

(tracks2 <- shifted |>
  filter(treatment %in% c('Bku', 'Bsu', 'Bst')) |>
  mutate(
    treatment = case_when(
      treatment == "Bku" ~ "*B. kuhniana*",
      treatment == "Bsu" ~ "*B. sudanica/<br>B. pfeifferi*",
      treatment == "Bst" ~ "*B. straminea*",
    ),
    treatment = factor(
      treatment,
      levels = c(
        "*B. sudanica/<br>B. pfeifferi*",
        "*B. straminea*",
        "*B. kuhniana*"
      )
    )
  ) |>
  ggplot(aes(x = x, y = y, color = frame / 8 / 60)) +
  geom_path(
    aes(group = interaction(experiment, subparticle)),
    linewidth = 0.25
  ) +
  # geom_point(size = 0.5) +
  facet_wrap(facets = vars(treatment), ncol = 3, scales = 'free_y') +
  scale_color_viridis_c(option = 'mako') +
  scale_x_continuous(
    breaks = seq(min(shifted$x), max(shifted$x), length.out = 5),
    expand = c(0, 0),
    labels = c("0", "20.25", "40.5", "60.75", "81")
  ) +
  scale_y_continuous(
    breaks = seq(min(shifted$y), max(shifted$y), length.out = 3),
    expand = c(0, 0),
    labels = c("", "", "")
  ) +
  labs(color = "Time\n(minutes)", x = ("X (mm)"), y = ("Y")) +
  theme_minimal() +
  theme(
    legend.position = "empty",
    strip.text.x = element_blank(),
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    strip.background = element_blank()
  ) +
  NULL)

save_plot(
  here("Fig3", "plots", "subplots", "Fig3D-bottom.pdf"),
  tracks2,
  base_width = 7.5,
  base_height = 3.5,
  device = cairo_pdf
)


########################################################
################### PLOT STACKED BAR ###################
########################################################

bins <- shifted |>
  ungroup() |>
  mutate(
    frame_bin = case_when(
      frame < 4800 * 2 ~ "0-20 min.",
      frame >= 4800 * 2 & frame < 4800 * 4 ~ "20-40 min.",
      frame >= 4800 * 4 ~ "40-60 min.",
    )
  ) |>
  # merge the APW treatment for plotting purposes
  mutate(
    treatment = case_when(
      str_detect(treatment, "APW") ~ "Control",
      treatment == "Bku" ~ "*B. kuhniana*",
      treatment == "Bga" ~ "*B. glabrata*",
      treatment == "Bsu" ~ "*B. sudanica/<br>B. pfeifferi*",
      treatment == "Bst" ~ "*B. straminea*",
    ),
    treatment = factor(
      treatment,
      levels = c(
        'Control',
        "*B. glabrata*",
        "*B. sudanica/<br>B. pfeifferi*",
        "*B. straminea*",
        "*B. kuhniana*"
      )
    )
  ) |>
  mutate(
    x_bin = case_when(
      x < -3000 ~ 'Start',
      x > -3000 & x < 3000 ~ "Middle",
      x > 3000 ~ 'Cue'
    ),
    x_bin = factor(
      x_bin,
      levels = c('Start', 'Middle', 'Cue')
    ),
    frame_bin = factor(
      frame_bin,
      levels = c('0-20 min.', '20-40 min.', "40-60 min.")
    )
  ) |>
  group_by(treatment, frame_bin, x_bin) %>%
  summarise(
    n = n()
  ) %>%
  mutate(pct = n / sum(n) * 100)

plot_stacked_bar <- function(df, cols = 2) {
  ggplot(df) +
    geom_col(aes(x = frame_bin, y = pct, fill = forcats::fct_rev(x_bin))) +
    coord_flip() +
    facet_wrap(~treatment, ncol = cols) +
    labs(x = 'Time bins', y = "Percent tracks in region", fill = "Region") +
    scale_fill_manual(
      values = rev(c("#FBDCE2", "#DA93C9", "#AE55BC")),
      guide = guide_legend(
        keywidth = 0.5,
        keyheight = 0.5
      )
    ) +
    scale_x_discrete(limits = rev) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_half_open() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      strip.text.x = element_markdown(size = 8),
      axis.text.x = element_markdown(size = 8, hjust = 1),
      axis.text.y = element_markdown(size = 8),
      axis.title = element_text(size = 9),
      strip.background = element_blank()
    )
}

(stacked_bar_1 <- bins %>%
  filter(treatment %in% c('Control', "*B. glabrata*")) |>
  plot_stacked_bar())

save_plot(
  here("Fig3", "plots", "subplots", "Fig3B-top.pdf"),
  stacked_bar_1,
  base_width = 7.5,
  base_height = 3.5,
  device = cairo_pdf
)

(stacked_bar_2 <- bins %>%
  filter(!treatment %in% c('Control', "*B. glabrata*")) |>
  plot_stacked_bar(cols = 3))

save_plot(
  here("Fig3", "plots", "subplots", "Fig3D-top.pdf"),
  stacked_bar_2,
  base_width = 7.5,
  base_height = 3.5,
  device = cairo_pdf
)

######################################################
################### CHUNK ANALYSIS ###################
######################################################

# nest_cols <- c("date", "experiment", "well_row", "well_col", "subparticle")
# nested <- quick_nest(chunked_data, nest_cols)

# subtrack_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 8,
#   chunk_size = 500,
#   pixel_to_mm = 126.5
# )

zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig3", "data"), files = 'skinny_subtrack_summary.rds')

subtrack_summary <- read_rds(here(
  "Fig3",
  "data",
  "skinny_subtrack_summary.rds"
))

subtrack_summary_region <- subtrack_summary |>
  mutate(
    response = ifelse(mean_x > 0, 'cue', 'control'),
    arena = str_extract(experiment, "a[0-9]{2}"),
    treatment = ifelse(
      well_row == 2,
      "Control",
      ifelse(
        arena == "a01" & well_row == 1,
        "Bku",
        ifelse(
          arena == "a01" & well_row == 0,
          "Bgl",
          ifelse(
            arena == "a02" & well_row == 1,
            "Bsu",
            ifelse(arena == "a02" & well_row == 0, "Bst", NA_character_)
          )
        )
      )
    ),
    .before = data,
  )

feature_cols <- subtrack_summary_region %>%
  ungroup() |>
  select(-(date:sd_y)) %>%
  select_if(is.numeric) %>%
  names()

source(here("utils", "model_utils.R"))

bgl_results <- map_dfr(
  feature_cols,
  ~ fit_model(.x, subtrack_summary_region |> filter(treatment == 'Bgl'))
) |>
  mutate(treatment = 'Bgl', .before = term)

bsu_results <- map_dfr(
  feature_cols,
  ~ fit_model(.x, subtrack_summary_region |> filter(treatment == 'Bsu'))
) |>
  mutate(treatment = 'Bsu', .before = term)

bku_results <- map_dfr(
  feature_cols,
  ~ fit_model(.x, subtrack_summary_region |> filter(treatment == 'Bku'))
) |>
  mutate(treatment = 'Bku', .before = term)

bst_results <- map_dfr(
  feature_cols,
  ~ fit_model(.x, subtrack_summary_region |> filter(treatment == 'Bst'))
) |>
  mutate(treatment = 'Bst', .before = term)

feature_track_summary <- subtrack_summary_region |>
  pivot_longer(
    cols = frame_start:curv_q90,
    names_to = "feature",
    values_to = "value"
  ) |>
  group_by(feature) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

results <- bind_rows(bgl_results, bsu_results, bku_results, bst_results) |>
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    sig = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ''
    )
  ) |>
  left_join(feature_track_summary |> select(feature, sd)) |>
  mutate(
    cohens_d = estimate / sd,
    se_cohens_d = std.error / sd
  ) |>
  arrange(desc(cohens_d))

feature_order <- results |>
  filter(str_detect(term, "cue")) |>
  group_by(feature) |>
  summarise(mean_cohens_d = mean(cohens_d, na.rm = TRUE), .groups = 'drop') |>
  arrange(mean_cohens_d) |>
  pull(feature)

feature_labels <- c(
  "directional_persistence" = "Directional persistence",
  "heading_autocorr" = "Heading autocorrelation",
  "angvel_q10" = "Angular velocity (10th percentile)",
  "accel_q10" = "Acceleration (10th percentile)",
  "angacc_q10" = "Angular acceleration (10th percentile)",
  "convex_hull_area" = "Convex hull area",
  "speed_autocorr" = "Speed autocorrelation",
  "net_displacement" = "Net displacement",
  "angvel_q50" = "Angular velocity (50th percentile)",
  "speed_q10" = "Speed (10th percentile)",
  "angacc_q50" = "Angular acceleration (50th percentile)",
  "movement_efficiency" = "Movement efficiency",
  "angular_acceleration_mean" = "Mean angular acceleration",
  "jerk_mean" = "Mean jerk",
  "radius_of_gyration" = "Radius of gyration",
  "acceleration_mean" = "Mean acceleration",
  "bbox_aspect_ratio" = "Bounding box aspect ratio",
  "angular_velocity_mean" = "Mean angular velocity",
  "turning_bias" = "Turning bias",
  "accel_q50" = "Acceleration (50th percentile)",
  "mean_heading" = "Mean heading",
  "speed_mean" = "Mean speed",
  "speed_q50" = "Speed (50th percentile)",
  "curv_q10" = "Curvature (10th percentile)",
  "path_length" = "Path length",
  "curv_q50" = "Curvature (50th percentile)",
  "fractal_dimension" = "Fractal dimension",
  "mean_curvature" = "Mean curvature",
  "curv_q90" = "Curvature (90th percentile)",
  "speed_q90" = "Speed (90th percentile)",
  "speed_var" = "Speed variance",
  "acceleration_var" = "Acceleration variance",
  "sinuosity" = "Sinuosity",
  "tortuosity" = "Tortuosity",
  "angular_velocity_var" = "Angular velocity variance",
  "angacc_q90" = "Angular acceleration (90th percentile)",
  "total_turn" = "Total turn",
  "angvel_q90" = "Angular velocity (90th percentile)",
  "accel_q90" = "Acceleration (90th percentile)",
  "heading_variance" = "Heading variance",
  "max_dist_from_start" = "Max distance from start",
  "speed_max" = "Max speed",
  "straightness" = "Straightness"
)

(results_plot <- results |>
  drop_na() |>
  mutate(
    feature = factor(feature, levels = names(feature_labels)),
    point_shape = ifelse(p_adj < 0.05, 16, NA)
  ) |>
  ggplot() +
  geom_rect(
    data = tibble(
      class = c('large', 'medium', 'small', 'medium', 'large'),
      xmin = c(-0.88, -0.5, -0.2, 0.2, 0.5),
      xmax = c(-0.5, -0.2, 0.2, 0.5, 0.6),
      ymin = -Inf,
      ymax = Inf
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = class),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_pointrange(aes(
    x = cohens_d,
    y = feature,
    color = treatment,
    shape = point_shape,
    xmin = cohens_d - se_cohens_d,
    xmax = cohens_d + se_cohens_d
  )) +
  scale_x_continuous(
    expand = expansion(mult = 0),
    breaks = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5)
  ) +
  scale_y_discrete(labels = feature_labels) +
  scale_color_manual(
    labels = c(
      "*B. glabrata*",
      "*B. kuhniana*",
      "*B. straminea*",
      "*B. sudanica/B. pfeifferi*"
    ),
    values = c("#EF476F", "#FFD166", "#06D6A0", "#118AB2")
  ) +
  scale_shape_identity() +
  scale_fill_manual(values = c('grey90', 'grey80', 'grey70')) +
  labs(x = "Standardized effect size", y = "Feature", color = "Cue source") +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 6),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_markdown(size = 8),
    legend.position = 'bottom'
  ) +
  NULL)

save_plot(
  here("Fig3", "plots", "subplots", "Fig3E.pdf"),
  results_plot,
  base_width = 10,
  base_height = 6,
  device = cairo_pdf
)


######################################################
################### FINAL PLOTTING ###################
######################################################

msa <- ggdraw() +
  draw_image(magick::image_read_pdf(here(
    "Fig3",
    "plots",
    "subplots",
    "Fig3C.pdf"
  )))

model <- ggdraw() +
  draw_image(magick::image_read_pdf(here(
    "Fig3",
    "plots",
    "subplots",
    "Fig3A.pdf"
  )))

library(patchwork)

(top_right <- (stacked_bar_1 + theme(legend.position = 'right')) /
  (tracks1 +
    theme_void() +
    theme(
      axis.title.y = element_blank(),
      strip.text.x = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
    )) +
  plot_layout(heights = c(2, 1.25)))

top_left <- plot_grid(model, msa, nrow = 2, labels = c("A", "C"), rel_heights = c(1, 1))

top <- plot_grid(
  top_left,
  top_right,
  nrow = 1,
  rel_widths = c(1, 1.25),
  labels = c("A", "B")
)

middleish <- (stacked_bar_2 + theme(legend.position = 'empty')) /
  (tracks2 +
    theme_void() +
    theme(axis.title.y = element_blank(), legend.position = 'empty', strip.text.x = element_blank())) +
  plot_layout(heights = c(2, 1))

final_plot <- plot_grid(
  top,
  middleish,
  results_plot,
  align = 'v',
  axis = 'r',
  labels = c("", "D", "E"),
  nrow = 3,
  rel_heights = c(1.2, 1, 2)
)

save_plot(
  here(
    "Fig3",
    "plots",
    "Fig3.pdf"
  ),
  final_plot,
  base_width = 7.5,
  base_height = 9,
  device = cairo_pdf
)

save_plot(
  here(
    "Fig3",
    "plots",
    "Fig3.png"
  ),
  final_plot,
  base_width = 7.5,
  base_height = 9,
  bg = 'white'
)
