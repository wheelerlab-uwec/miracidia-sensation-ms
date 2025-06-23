library(tidyverse)
library(here)
library(tidymodels)
library(ggtext)
library(cowplot)

usethis::use_air(vscode = TRUE)

# functions ---------------------------------------------------------------

source(here("utils", "helper_functions.R"))

# import and tidy ---------------------------------------------------------

# 4.75 GB file of raw tracking/coordinate data
# NOTE: can skip this and read the summary in line 54 or uncomment and run to see raw tracks
# zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig1", "data"), files = 'double_agar.rds')

# double_agar <- read_rds(here(
#   "Fig1",
#   "data",
#   "double_agar.rds"
# ))

# # keep any track > 40 frames (8 FPS means ~5 seconds)
# frame_filter <- double_agar |>
#   group_by(file, date, video, experiment, particle) |>
#   tally() |>
#   filter(n > 40)

# double_agar_filtered <- frame_filter |>
#   select(-n) |>
#   left_join(double_agar)

# # split all tracks into 5 sec chunks
# chunked_data <- double_agar_filtered |>
#   group_by(file, date, video, experiment, particle) |>
#   arrange(frame, .by_group = TRUE) |>
#   group_split() |>
#   map_dfr(~ split_trajectory(.x, frame_rate = 8, chunk_duration_sec = 5))

# nest_cols <- c("file", "date", "video", "experiment", "particle")

# nested <- quick_nest(chunked_data, nest_cols)

# # calculate features for each track (particle == track in trackpy terminology)
# track_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 8,
#   chunk_size = 500
# )

# NOTE: start here
zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig1", "data"), files = 'double_track_summary.rds')

track_summary <- read_rds(here("Fig1", "data", 'double_track_summary.rds'))

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


#####################################################
################### PLOT HEATMAPS ###################
#####################################################

plot_heatmap <- function(df, fill_limits = c(0, 1)) {
  ggplot(df) +
    geom_hex(
      aes(x, y, fill = after_stat(density)),
      binwidth = c(350, 350),
      stat = "binhex"
    ) +
    annotate(
      "point",
      x = -192.5,
      y = 1712,
      size = 5,
      fill = "orange",
      color = "white",
      alpha = 0.50,
      shape = 21
    ) +
    # these hard-coded values were derived from overlays of the diffusion experiments
    # attempts were made to make the responses the same size and equidistant from the entry pore
    # control response
    annotate(
      "rect",
      xmin = min(shifted$x) + 25,
      xmax = min(shifted$x) + 3693,
      ymin = 0,
      ymax = 3400,
      color = "steelblue",
      fill = NA,
      linetype = "dashed",
      linewidth = 2
    ) +
    # scw response
    annotate(
      "rect",
      xmin = min(shifted$x) + 7129 - 400,
      xmax = max(shifted$x),
      ymin = 0,
      ymax = 3400,
      color = "springgreen4",
      fill = NA,
      linetype = "dashed",
      linewidth = 2
    ) +
    # scw agar
    annotate(
      "rect",
      xmin = max(shifted$x) - 1180 - 1000,
      xmax = max(shifted$x) - 1180,
      ymin = 1200,
      ymax = 1200 + 1000,
      fill = "springgreen4",
      alpha = .50,
      color = "white"
    ) +
    # control agar
    annotate(
      "rect",
      xmin = min(shifted$x) + 1264,
      xmax = min(shifted$x) + 1264 + 1000,
      ymin = 1200,
      ymax = 1200 + 1000,
      fill = "steelblue",
      alpha = .50,
      color = "white"
    ) +
    annotate(
      "richtext",
      x = 3200,
      y = 1700,
      label = "***Snail<br>cue***",
      color = "white",
      size = 2,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    annotate(
      "richtext",
      x = -3700,
      y = 1700,
      label = "***Control***",
      color = "white",
      size = 2,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    coord_equal() +
    scale_x_continuous(
      breaks = seq(min(shifted$x), max(shifted$x), length.out = 5),
      expand = c(0, 0),
      labels = c("0", "20.25", "40.5", "60.75", "81")
    ) +
    scale_y_continuous(
      breaks = seq(min(shifted$y), max(shifted$y), length.out = 5),
      expand = c(0, 0),
      labels = c("", "6.75", "13.5", "20.25", "27")
    ) +
    scale_fill_viridis_c(option = "magma", limits = fill_limits) +
    labs(color = "Time\n(minutes)", x = ("X (mm)"), y = ("Y (mm)")) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = 'black'),
      legend.position = "empty",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
    ) +
    NULL
}

# extract computed densities from a combined plot so separate plots use the same color scale
combined_plot <- ggplot(shifted, aes(x, y)) +
  geom_hex(binwidth = c(350, 350), stat = "binhex") +
  facet_grid(scw_type ~ diffusion)

density_vals <- ggplot_build(combined_plot)$data[[1]]$density
max_density <- max(density_vals, na.rm = TRUE)
fill_limits <- c(0, max_density)

(heatmap_1 <- shifted |>
  filter(diffusion == "1 hr.") |>
  plot_heatmap(fill_limits = fill_limits))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1E-left.pdf"),
  heatmap_1,
  base_width = 6,
  base_height = 3
)

(heatmap_2 <- shifted |>
  filter(diffusion == "2 hr." & scw_type == 'lyo') |>
  plot_heatmap(fill_limits = fill_limits))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1D-left.pdf"),
  heatmap_2,
  base_width = 6,
  base_height = 3
)

(heatmap_3 <- shifted |>
  filter(diffusion == "2 hr." & scw_type == 'raw') |>
  plot_heatmap(fill_limits = fill_limits))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1C-left.pdf"),
  heatmap_3,
  base_width = 6,
  base_height = 3
)


######################################################
################### PLOT DENSITIES ###################
######################################################

plot_density <- function(df, diff_rects) {
  ggplot(df) +
    # control response
    annotate(
      "rect",
      xmin = min(shifted$x) + 25,
      xmax = min(shifted$x) + 3693,
      ymin = 0,
      ymax = Inf,
      fill = "steelblue",
      alpha = .5
    ) +
    # scw response
    annotate(
      "rect",
      xmin = min(shifted$x) + 7129 - 400,
      xmax = max(shifted$x),
      ymin = 0,
      ymax = Inf,
      fill = "springgreen4",
      alpha = .5
    ) +
    annotate(
      "richtext",
      x = 3250,
      y = 1.5e-4,
      label = "***Snail cue***",
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt"),
      color = "white",
      size = 3
    ) +
    annotate(
      "richtext",
      x = -3750,
      y = 1.5e-4,
      label = "***Control***",
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(rep(0, 4), "pt"),
      color = "white",
      size = 3
    ) +
    geom_density(aes(x = x), adjust = 2, linewidth = 2) +
    geom_density(
      aes(x = x, color = video),
      adjust = 4,
      linewidth = 1,
      alpha = 0.25
    ) +
    ggarrow::geom_arrow(
      data = diff_rects,
      aes(x = x, y = y),
      color = 'indianred',
    ) +
    scale_color_manual(
      values = alpha(
        c("grey1", "grey20", "grey40", "grey60", "grey80", "grey90"),
        0.4
      )
    ) +
    scale_x_continuous(
      limits = c(min(shifted$x), max(shifted$x)),
      breaks = seq(min(shifted$x), max(shifted$x), length.out = 5),
      expand = c(0, 0),
      labels = c("0", "20.25", "40.5", "60.75", "81")
    ) +
    scale_y_continuous(
      labels = scales::scientific,
      n.breaks = 4
    ) +
    # facet_grid(rows = vars(diffusion)) +
    labs(x = "X (mm)", y = "Density") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.ticks = element_line(),
      legend.position = "none"
    ) +
    NULL
}

diff_rects_1 <- tibble(
  diffusion = c("1 hr.", "1 hr."),
  x = c(min(shifted$x) + 3693, min(shifted$x) + 7129 - 400),
  y = c(0.925e-04, 1.2e-04),
)

diff_rects_2 <- tibble(
  diffusion = c("2 hr.", "2 hr."),
  x = c(min(shifted$x) + 3693, min(shifted$x) + 7129 - 400),
  y = c(1.1e-04, 1.55e-04),
)

diff_rects_3 <- tibble(
  diffusion = c("2 hr.", "2 hr."),
  x = c(min(shifted$x) + 3693, min(shifted$x) + 7129 - 400),
  y = c(1.e-04, 1.42e-04),
)

(density_1 <- shifted |>
  filter(diffusion == '1 hr.') |>
  plot_density(diff_rects = diff_rects_1))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1E-right.pdf"),
  density_1,
  base_width = 6,
  base_height = 3
)

(density_2 <- shifted |>
  filter(diffusion == "2 hr." & scw_type == 'lyo') |>
  plot_density(diff_rects = diff_rects_2))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1D-right.pdf"),
  density_2,
  base_width = 6,
  base_height = 3
)

(density_3 <- shifted |>
  filter(diffusion == "2 hr." & scw_type == 'raw') |>
  plot_density(diff_rects = diff_rects_3))

save_plot(
  here("Fig1", "plots", "subplots", "Fig1C-right.pdf"),
  density_3,
  base_width = 6,
  base_height = 3
)


#########################################################
################### SUBTRACK ANALYSIS ###################
#########################################################

# NOTE: can skip this and read the summary in line 416
nest_cols <- c("file", "date", "video", "experiment", "particle", "subparticle")

chunked_nested <- quick_nest(
  shifted |> select(-(frame_start:curv_q90)),
  nest_cols
)

# calculate features for the track chunks
subtrack_summary <- calculate_track_features_parallel(
  chunked_nested,
  fps = 8,
  chunk_size = 500,
  pixel_to_mm = 126.5
)

zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig1", "data"), files = 'double_subtrack_summary.rds')

subtrack_summary <- read_rds(here("Fig1", "data", 'double_subtrack_summary.rds')) |>
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
    ),
    # response (generic variable so we can use fit_model across experiments) == region of interest
    response = case_when(
      mean_x > 1247.013 ~ 'cue',
      mean_x < -1788.987 ~ 'control',
      TRUE ~ 'middle'
    ),
    .before = data
  ) |>
  # don't include the middle region in the modeling
  filter(response != 'middle')

source(here("utils", "model_utils.R"))

feature_cols <- subtrack_summary %>%
  ungroup() |>
  select(-(file:sd_y)) %>%
  select_if(is.numeric) %>%
  names()

d1_results <- map_dfr(
  feature_cols,
  ~ fit_model(.x, subtrack_summary |> filter(diffusion == '1 hr.'))
) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

d2l_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    subtrack_summary |> filter(diffusion == "2 hr." & scw_type == 'lyo')
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

d2r_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    subtrack_summary |> filter(diffusion == "2 hr." & scw_type == 'raw')
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

feature_summary <- subtrack_summary |>
  pivot_longer(
    cols = frame_start:curv_q90,
    names_to = "feature",
    values_to = "value"
  ) |>
  group_by(scw_type, diffusion, feature) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

results <- bind_rows(
  d1_results |> mutate(diffusion = '1 hr.', scw_type = 'lyo'),
  d2l_results |> mutate(diffusion = '2 hr.', scw_type = 'lyo'),
  d2r_results |> mutate(diffusion = '2 hr.', scw_type = 'raw')
) |>
  mutate(
    sig = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ''
    )
  ) |>
  left_join(feature_summary |> select(scw_type, diffusion, feature, sd)) |>
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
  filter(str_detect(term, "cue")) |>
  drop_na() |>
  mutate(
    feature = factor(feature, levels = feature_order),
    linestart = cohens_d - se_cohens_d,
    lineend = cohens_d + se_cohens_d,
    cohens_d = ifelse(p_adj > 0.05, NA, cohens_d)
  ) |>
  ggplot() +
  geom_rect(
    data = tibble(
      class = c('medium', 'small', 'medium'),
      xmin = c(-0.45, -0.2, 0.2),
      xmax = c(-0.2, 0.2, 0.45),
      ymin = -Inf,
      ymax = Inf
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = class),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_segment(aes(
    x = linestart,
    xend = lineend,
    y = feature,
    color = diffusion,
    linetype = scw_type,
  )) +
  geom_point(
    aes(
      x = cohens_d,
      y = feature,
      color = diffusion,
      shape = scw_type
    ),
    size = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0),
    breaks = c(-0.4, -0.2, 0, 0.2, 0.4)
  ) +
  scale_y_discrete(labels = feature_labels) +
  scale_color_manual(
    labels = c("None", "1 hr."),
    values = c("indianred", "steelblue")
  ) +
  scale_linetype_manual(
    values = c("solid", "dashed"),
    labels = c("Lyophilized", "Raw")
  ) +
  scale_shape_manual(
    values = c(16, 17),
    labels = c("Lyophilized", "Raw"),
  ) +
  scale_fill_manual(values = c('grey90', 'grey80')) +
  labs(
    x = "Standardized effect size",
    y = "Feature",
    color = "Pre-diffusion",
    shape = "SCW type",
    linetype = "SCW type"
  ) +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_markdown(size = 8),
    legend.position = 'right'
  ) +
  NULL)

save_plot(
  here("Fig1", "plots", "subplots", "Fig1F.pdf"),
  results_plot,
  base_width = 10,
  base_height = 4
)


######################################################
################### FINAL PLOTTING ###################
######################################################

model <- ggdraw() +
  draw_image(magick::image_read_pdf(here(
    "Fig1",
    "plots",
    "subplots",
    "Fig1A.pdf"
  )))

diffusion <- ggdraw() +
  draw_image(magick::image_read_pdf(here("Fig1", "plots", "subplots", "Fig1B.pdf")))

top <- plot_grid(
  model,
  diffusion,
  nrow = 1,
  labels = c("A", "B"),
  rel_widths = c(0.75, 1)
)

middle <- plot_grid(
  heatmap_3 +
    ggtitle("1 hr. pre-diffusion (raw)") +
    theme(plot.title = element_text(hjust = 0.5)),
  density_3,

  heatmap_2 +
    ggtitle("1 hr. pre-diffusion (lyophilized)") +
    theme(plot.title = element_text(hjust = 0.5)),
  density_2,

  heatmap_1 +
    ggtitle("No pre-diffusion (lyophilized)") +
    theme(plot.title = element_text(hjust = 0.5)),
  density_1,

  nrow = 3,
  ncol = 2,
  rel_widths = c(0.7, 1),
  align = "hv",
  axis = "tblr",
  labels = c("C", "", "D", "", "E")
)

final_plot <- plot_grid(
  top,
  middle,
  results_plot,
  labels = c("", "", "F"),
  nrow = 3,
  rel_heights = c(0.45, 1.25, 1)
)

save_plot(
  here("Fig1", "plots", "Fig1.pdf"),
  final_plot,
  base_width = 10,
  base_height = 12,
  device = cairo_pdf
)

save_plot(
  here("Fig1", "plots", "Fig1.png"),
  final_plot,
  base_width = 10,
  base_height = 12,
  bg = 'white'
)
