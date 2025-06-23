library(tidyverse)
library(cowplot)
library(here)
library(ggtext)

# functions ---------------------------------------------------------------

source(here("../../invision/double_agar/helper_functions.R"))

# import ------------------------------------------------------------------

# frame_filter <- tracks |>
#   group_by(file, experiment, particle) |>
#   tally() |>
#   # keep any track > 5 seconds (5 x 15 fps = 75 frames)
#   filter(n > 75)

# filtered <- frame_filter |>
#   select(-n) |>
#   left_join(tracks) %>%
#   select(-video)

# chunked_data <- filtered |>
#   group_by(file, date, experiment, particle) |>
#   arrange(frame, .by_group = TRUE) |>
#   group_split() |>
#   map_dfr(~ split_trajectory(.x, frame_rate = 15, chunk_duration_sec = 5))

# nest_cols <- c(
#   "date",
#   "file",
#   "experiment",
#   "timepoint",
#   "side",
#   "particle",
#   "subparticle"
# )

# nested <- quick_nest(chunked_data, nest_cols)

# screening_subtrack_summary <- calculate_track_features_parallel(
#   nested,
#   fps = 15,
#   chunk_size = 1000
# )

# NOTE: start here
zen4R::download_zenodo("10.5281/zenodo.15713500", path = here("Fig4", "data"), files = 'screening_subtrack_summary.rds')
screening_subtrack_summary.rds <- readRDS(here("Fig4", "data", "screening_subtrack_summary.rds"))

shifted <- screening_subtrack_summary.rds |>
  unnest(c(data)) |>
  group_by(date, experiment, particle, subparticle) |>
  arrange(frame, .by_group = TRUE) |>
  group_by(experiment) |>
  mutate(
    x = case_when(
      side == "left" ~ -x, # flip left video horizontally and place on negative x-axis
      side == "right" ~ max(x) - x # flip right video horizontally but keep on positive side
    )
  ) |>
  mutate(experiment = str_sub(experiment, end = -10))

split_wells <- shifted |>
  mutate(
    well_row = case_when(
      well_row == 0 ~ "A",
      well_row == 1 ~ "B",
      well_row == 2 ~ "C",
      well_row == 3 ~ "D",
      well_row == 4 ~ "E",
      well_row == 5 ~ "F",
    ),
    well_col = case_when(
      well_col == 0 & side == "left" ~ "01",
      well_col == 1 & side == "left" ~ "02",
      well_col == 2 & side == "left" ~ "03",
      well_col == 3 & side == "left" ~ "04",
      well_col == 4 & side == "left" ~ "05",
      well_col == 5 & side == "left" ~ "06",
      well_col == 6 & side == "left" ~ "07",
      well_col == 7 & side == "left" ~ "08",
      well_col == 0 & side == "right" ~ "09",
      well_col == 1 & side == "right" ~ "10",
      well_col == 2 & side == "right" ~ "11",
      well_col == 3 & side == "right" ~ "12",
      well_col == 4 & side == "right" ~ "13"
    ),
    well = paste0(well_row, well_col),
    treatment = case_when(
      well_col == "01" ~ "APW",
      well_col == "02" ~ "5 mM MgCl<sub>2</sub>",
      well_col == "03" ~ "10 mM MgCl<sub>2</sub>",
      well_col == "04" ~ "0.2 ppm *B. glabrata* P12",
      well_col == "05" ~ "0.4 ppm *B. glabrata* P12",
      well_col == "06" ~ "0.2 ppm *B. kuhniana* P12",
      well_col == "07" ~ "0.4 ppm *B. kuhniana* P12",
      well_col == "08" ~ "0.2 ppm *B. sudanica* P12",
      well_col == "09" ~ "0.4 ppm *B. sudanica* P12",
      well_col == "10" ~ "0.2 ppm *B. straminea* P12",
      well_col == "11" ~ "0.4 ppm *B. straminea* P12",
      well_col == "12" ~ "1X lSCW",
      well_col == "13" ~ "2X lSCW",
    ),
    treatment = factor(
      treatment,
      levels = c(
        "APW",
        "5 mM MgCl<sub>2</sub>",
        "10 mM MgCl<sub>2</sub>",
        "1X lSCW",
        "2X lSCW",
        "0.2 ppm *B. glabrata* P12",
        "0.4 ppm *B. glabrata* P12",
        "0.2 ppm *B. kuhniana* P12",
        "0.4 ppm *B. kuhniana* P12",
        "0.2 ppm *B. sudanica* P12",
        "0.4 ppm *B. sudanica* P12",
        "0.2 ppm *B. straminea* P12",
        "0.4 ppm *B. straminea* P12"
      )
    )
  )

track_prep <- split_wells %>%
  filter(
    timepoint %in% c("30 min", "60 min"),
    well_row != "A",
  ) %>%
  mutate(timepoint = as.numeric(str_remove(timepoint, ' min'))) %>%
  # move over a bit
  mutate(
    x = case_when(
      side == 'right' ~ x - 175,
      TRUE ~ x
    ),
    y = case_when(
      side == 'right' ~ y + 115,
      TRUE ~ y
    )
  )

ex_plot <- track_prep |>
  filter(
    well %in% c("B01", "B02", "B03"),
    date == '20250410',
    timepoint == 30
  ) |>
  ggplot() +
  geom_path(
    aes(
      x = x,
      y = y,
      group = interaction(well, particle),
      color = frame / 15 / 60
    ),
    linewidth = 0.25
  ) +
  facet_wrap(vars(timepoint), nrow = 2) +
  scale_color_viridis_c(
    option = "mako",
    guide = guide_colourbar(barheight = 3, barwidth = 0.5)
  ) +
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
  coord_equal() +
  labs(x = "X position (mm)", y = "Y position (mm)", color = "Minute") +
  theme_void() +
  theme(
    # aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    # legend.text = element_text(angle = 45, hjust = 1),
    strip.text = element_blank()
  )
ex_plot

# subtrack analysis ----------------------------------------------------

source(here("Fig4", "model_utils.R"))

analysis <- split_wells |>
  select(-x, -y, -frame) |>
  distinct() |>
  filter(
    sd_x & sd_y > 1
  )

feature_cols <- analysis %>%
  ungroup() |>
  select(-(file:sd_y)) %>%
  select_if(is.numeric) %>%
  names()

t0_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    analysis |> filter(timepoint == "0 min"),
    control_level = "APW"
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "holm"))

t30_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    analysis |> filter(timepoint == "30 min"),
    control_level = "APW"
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "holm"))

t60_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    analysis |> filter(timepoint == "60 min"),
    control_level = "APW"
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "holm"))

t90_results <- map_dfr(
  feature_cols,
  ~ fit_model(
    .x,
    analysis |> filter(timepoint == "90 min"),
    control_level = "APW"
  )
) |>
  mutate(p_adj = p.adjust(p.value, method = "holm"))

feature_summary <- analysis |>
  pivot_longer(
    cols = frame_start:curv_q90,
    names_to = "feature",
    values_to = "value"
  ) |>
  group_by(timepoint, feature) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) |>
  drop_na()

results <- bind_rows(
  t0_results |> mutate(timepoint = "0 min"),
  t30_results |> mutate(timepoint = "30 min"),
  t60_results |> mutate(timepoint = "60 min"),
  t90_results |> mutate(timepoint = "90 min"),
) |>
  mutate(
    sig = case_when(
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      p_adj > 0.05 ~ ""
    ),
    treatment = str_remove(term, "treatment"),
  ) |>
  select(-term) |>
  left_join(feature_summary |> select(timepoint, feature, sd)) |>
  mutate(
    cohens_d = estimate / sd,
    se_cohens_d = std.error / sd
  ) |>
  arrange(desc(cohens_d))

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

(feature_plot <- results |>
  drop_na() |>
  filter(timepoint %in% c("30 min", "60 min")) |>
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "5 mM MgCl<sub>2</sub>",
        "10 mM MgCl<sub>2</sub>",
        "1X lSCW",
        "2X lSCW",
        "0.2 ppm *B. glabrata* P12",
        "0.4 ppm *B. glabrata* P12",
        "0.2 ppm *B. kuhniana* P12",
        "0.4 ppm *B. kuhniana* P12",
        "0.2 ppm *B. sudanica* P12",
        "0.4 ppm *B. sudanica* P12",
        "0.2 ppm *B. straminea* P12",
        "0.4 ppm *B. straminea* P12"
      )
    )
  ) |>
  ggplot() +
  geom_tile(aes(x = treatment, y = feature, fill = cohens_d)) +
  geom_text(
    aes(x = treatment, y = feature, label = sig),
    color = "grey10",
    size = 1.5,
    vjust = 0.5
  ) +
  facet_wrap(vars(timepoint), nrow = 1) +
  scale_y_discrete(labels = feature_labels) +
  scale_fill_distiller(
    limits = c(-1, 1),
    palette = "BrBG"
  ) +
  labs(fill = "Standardized<br>effect size", y = "Feature", x = "Treatment") +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    strip.text = element_markdown(size = 9),
    legend.title = element_markdown(size = 9),
    legend.text = element_text(size = 8),
    legend.position = 'bottom'
  ) +
  NULL)

supp_features <- results |>
  drop_na() |>
  filter(timepoint %in% c("0 min", "90 min")) |>
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "5 mM MgCl<sub>2</sub>",
        "10 mM MgCl<sub>2</sub>",
        "1X lSCW",
        "2X lSCW",
        "0.2 ppm *B. glabrata* P12",
        "0.4 ppm *B. glabrata* P12",
        "0.2 ppm *B. kuhniana* P12",
        "0.4 ppm *B. kuhniana* P12",
        "0.2 ppm *B. sudanica* P12",
        "0.4 ppm *B. sudanica* P12",
        "0.2 ppm *B. straminea* P12",
        "0.4 ppm *B. straminea* P12"
      )
    )
  ) |>
  ggplot() +
  geom_tile(aes(x = treatment, y = feature, fill = cohens_d)) +
  geom_text(
    aes(x = treatment, y = feature, label = sig),
    color = "grey10",
    size = 1.5,
    vjust = 0.5
  ) +
  facet_wrap(vars(timepoint), nrow = 1) +
  scale_y_discrete(labels = feature_labels) +
  scale_fill_distiller(
    limits = c(-1, 1),
    palette = "BrBG"
  ) +
  labs(fill = "Standardized<br>effect size", y = "Feature", x = "Treatment") +
  theme_half_open() +
  theme(
    axis.text.x = element_markdown(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_markdown(size = 8),
    axis.title = element_text(size = 9),
    strip.text = element_markdown(size = 9),
    legend.title = element_markdown(size = 9),
    legend.text = element_text(size = 8),
    legend.position = 'bottom'
  ) +
  NULL


save_plot(
  here("Fig4", "plots", "S3_Fig.pdf"),
  supp_features,
  base_width = 6,
  base_height = 6
)

save_plot(
  here("Fig4", "plots", "S3_Fig.png"),
  supp_features,
  base_width = 7.5,
  base_height = 7.5,
  bg = 'white'
)

# final plotting

model <- ggdraw() +
  draw_image(magick::image_read_pdf(here(
    "Fig4",
    "plots",
    "subplots",
    "Fig4A.pdf"
  )))

top <- plot_grid(model, ex_plot, ncol = 2, labels = c('A', 'B'), rel_widths = c(1.25, 1.5))

merge <- plot_grid(
  top,
  feature_plot + theme(legend.position = 'right'),
  nrow = 2,
  rel_heights = c(0.4, 1.5),
  labels = c('', 'C')
)

save_plot(
  here("Fig4", "plots", "Fig4.pdf"),
  merge,
  base_width = 7.5,
  base_height = 7,
  device = cairo_pdf
)

save_plot(
  here("Fig4", "plots", "Fig4.png"),
  merge,
  base_width = 7.5,
  base_height = 7,
  bg = 'white'
)
