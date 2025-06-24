library(tidyverse)
library(here)
library(tidymodels)
library(ggtext)
library(cowplot)

usethis::use_air(vscode = TRUE)

# functions ---------------------------------------------------------------

source(here("utils", "helper_functions.R"))

# import and tidy ---------------------------------------------------------

# tracks with features extracted (see Fig1.R for details)
track_summary <- read_rds(here("Fig1", "data", "double_track_summary.rds"))

shifted <- track_summary |>
  unnest(c(data)) |>
  group_by(date, video, experiment, particle) |>
  arrange(frame, .by_group = TRUE) |>
  group_by(date, video, experiment) |>
  mutate(
    diffusion = case_when(
      date %in% c("20240111", "20240328", "20240627", "20240718", '20240725') ~ "2 hr.",
      date %in% c("20240808", "20240822", '20240829') ~ "1 hr."
    ),
    scw_type = case_when(
      date %in% c("20240111", "20240328") ~ 'raw',
      date %in%
        c(
          "20240627",
          "20240718",
          '20240725',
          "20240808",
          "20240822",
          '20240829'
        ) ~
        "lyo"
    )
  ) |>
  # filter out non-moving particles (false-positives)
  filter(
    speed_mean > 5,
    net_displacement > 50,
    sd_x > 5 & sd_y > 5
  )

apw_border <- min(shifted$x) + 3693
scw_border <- min(shifted$x) + 7129 - 400

border_crossers <- shifted |>
  filter(diffusion == '2 hr.') |>
  group_by(file, date, video, experiment, subparticle) |>
  slice(1, n()) |>
  select(file:particle, subparticle, x) |>
  mutate(position = rep(c("start", "stop"))) |>
  pivot_wider(names_from = position, values_from = x) |>
  mutate(
    category = case_when(
      start > apw_border & stop < apw_border ~ "out_in_APW",
      start < apw_border & stop > apw_border ~ "in_out_APW",
      start < scw_border & stop > scw_border ~ "out_in_SCW",
      start > scw_border & stop < scw_border ~ "in_out_SCW",
      TRUE ~ NA_character_
    )
  ) |>
  drop_na(category) |>
  select(-start, -stop) |>
  left_join(select(shifted, file:x, subparticle, frame))

border_crossers_features <- border_crossers |>
  rename(response = category) |>
  group_by(file, date, video, experiment, particle, subparticle, response) |>
  group_nest() |>
  calculate_track_features_parallel(
    fps = 8,
    pixel_to_mm = 126.5,
    chunk_size = 50
  )

border_crossers_features_summary <- border_crossers_features |>
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

feature_cols <- border_crossers_features %>%
  ungroup() |>
  select(-(file:sd_y)) %>%
  select_if(is.numeric) %>%
  names()

source(here("utils", "model_utils.R"))

out_in_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    border_crossers_features |> filter(str_detect(response, "^out_in"))
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

in_out_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    border_crossers_features |> filter(str_detect(response, "^in_out"))
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

border_crossers_results <- bind_rows(out_in_results, in_out_results) |>
  mutate(
    sig = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ''
    )
  ) |>
  left_join(border_crossers_features_summary |> select(feature, sd)) |>
  mutate(
    cohens_d = estimate / sd,
    se_cohens_d = std.error / sd
  ) |>
  arrange(desc(cohens_d)) |>
  drop_na(p_adj)

feature_order <- border_crossers_results |>
  drop_na(p_adj) |>
  group_by(feature) |>
  summarise(mean_cohens_d = mean(cohens_d, na.rm = TRUE), .groups = 'drop') |>
  arrange(mean_cohens_d) |>
  drop_na() |>
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

(border_crossers_cleveland <- border_crossers_results |>
  drop_na(cohens_d) |>
  mutate(
    point_shape = ifelse(p_adj <= 0.05, "sig", NA),
    feature = factor(feature, levels = feature_order),
  ) |>
  ggplot() +
  geom_rect(
    data = tibble(
      class = c('large', 'medium', 'small', 'medium', 'large'),
      xmin = c(-1, -0.8, -0.5, 0.5, 0.8),
      xmax = c(-0.8, -0.5, 0.5, 0.8, 1),
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
    color = term,
    xmin = cohens_d - se_cohens_d,
    xmax = cohens_d + se_cohens_d,
    shape = point_shape
  )) +
  scale_x_continuous(
    limits = c(-1, 1),
    expand = expansion(mult = c(0, 0)),
    breaks = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8)
  ) +
  scale_y_discrete(labels = feature_labels) +
  scale_color_manual(
    values = c("darkorange", "purple"),
    labels = c("In to out", "Out to in")
  ) +
  scale_fill_manual(values = c('grey90', 'grey80', 'grey70')) +
  labs(
    x = "Standardized effect size",
    y = "Feature",
    color = "Bording crossing direction"
  ) +
  guides(
    shape = "none"
  ) +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_markdown(size = 8),
    legend.position = 'bottom'
  ) +
  NULL)

save_plot(
  here("Fig2", "plots", "subplots", 'Fig2A.pdf'),
  border_crossers_cleveland,
  base_width = 6,
  base_height = 5
)

(border_crossers_plot <- border_crossers_features |>
  filter(str_detect(response, "SCW")) |>
  mutate(response = factor(response, levels = c("out_in_SCW", "in_out_SCW"))) |>
  unnest(data) |>
  ggplot() +
  geom_path(
    aes(x, y, color = angular_velocity_var, group = subparticle),
    arrow = arrow(type = "closed", length = unit(0.05, "inches")),
    linewidth = 0.5,
    alpha = 0.75
  ) +
  geom_vline(xintercept = scw_border, linetype = "dashed", color = "white") +
  scale_color_viridis_c(option = "mako") +
  facet_wrap(
    vars(response),
    scale = 'free',
    labeller = as_labeller(
      c(
        "out_in_SCW" = "Out to in - SCW",
        "in_out_SCW" = "In to out - SCW"
      )
    )
  ) +
  labs(x = ("X (px)"), y = "Y (px)", color = "Angular velocity\nvariance") +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = 'right'
  ) +
  NULL)

save_plot(
  here("Fig2", "plots", "subplots", 'Fig2B.pdf'),
  border_crossers_plot,
  base_width = 6,
  base_height = 5
)

save_plot(
  here("Fig2", "plots", 'Fig2.pdf'),
  plot_grid(
    border_crossers_cleveland,
    border_crossers_plot,
    nrow = 2,
    labels = c("A", "B"),
    rel_heights = c(1, 0.75)
  ),
  base_width = 6,
  base_height = 8
)

save_plot(
  here("Fig2", "plots", 'Fig2.png'),
  plot_grid(
    border_crossers_cleveland,
    border_crossers_plot,
    nrow = 2,
    labels = c("A", "B"),
    rel_heights = c(1, 0.75)
  ),
  base_width = 6,
  base_height = 8,
  bg = 'white'
)
