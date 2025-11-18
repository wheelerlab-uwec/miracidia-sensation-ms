library(here)
library(tidyverse)
library(multcomp)
library(cowplot)
library(ggtext)

df <- read_csv(here("Fig5", "data", "penetration_data.csv")) %>%
  mutate(
    infections = before - after,
    percent = infections / before,
    treatment = case_when(
      treatment == 'water' ~ 'Water',
      TRUE ~ treatment
    ),
    treatment = factor(treatment, levels = c('Water', "Bgla", "Bsud", "Bstr", "Bkun"))
  ) %>%
  mutate() %>%
  #filtering out preliminary data - didn't actually have pre-infection counts
  filter(experiment_date != '20250213', before != 10, )


summary <- df %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(percent),
    sd = sd(percent)
  )

comparisons <- list(c('Water', 'Bgla'), c('Water', 'Bsud'), c('Water', 'Bstr'), c('Water', 'Bkun'))

df$treatment <- relevel(df$treatment, ref = 'Water')

model <- lm(percent ~ treatment, data = df)
dt <- glht(model, linfct = mcp(treatment = "Dunnett")) |>
  broom::tidy() |>
  mutate(
    treatment = str_remove(contrast, ' - Water'),
    sig = case_when(
      adj.p.value < 0.0001 ~ '****',
      adj.p.value < 0.001 ~ '***',
      adj.p.value < 0.01 ~ '**',
      adj.p.value < 0.05 ~ '*',
      TRUE ~ 'ns'
    )
  )

plot_grid(
  ggpubr::ggqqplot(residuals(model)) + ggtitle("QQ Plot of Residuals")
)

labels <- c(
  Water = "Water",
  Bgla = "<img src='https://i.postimg.cc/B6ZpnzxZ/Screenshot-2025-06-05-at-3-39-25-PM.png'
    width='150' /><br>*B. glabrata*",
  Bsud = "<img src='https://i.postimg.cc/Px7FYNRN/Screenshot-2025-04-03-at-8-43-15-PM.png'
    width='150' /><br>*B. sudanica*<br>*B. pfeifferi*",
  Bstr = "<img src='https://i.postimg.cc/FsNnctPy/Screenshot-2025-04-03-at-8-44-31-PM.png'
    width='150' /><br>*B. straminea*",
  Bkun = "<img src='https://i.postimg.cc/jjH7y13f/Screenshot-2025-04-03-at-8-45-28-PM.png'
    width='150' /><br>*B. kuhniana*"
)


(plot <- df %>%
  ggplot(aes(x = treatment, y = percent)) +
  # geom_bar(stat = "summary", fun = "mean") +
  geom_boxplot(
    outlier.shape = NA,
    fill = 'lightgrey',
    color = 'grey30',
    alpha = 0.7,
    width = 0.6
  ) +
  ggbeeswarm::geom_quasirandom(size = 2, shape = 21, color = 'white', fill = 'steelblue') +
  # stat_summary(fun.data = "mean_cl_normal", color = 'indianred') +
  geom_text(
    data = dt,
    aes(
      x = treatment,
      y = 1.05,
      label = sig
    ),
    inherit.aes = FALSE,
    size = 4,
    color = 'black'
  ) +
  scale_x_discrete(
    name = NULL,
    labels = labels
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    expand = expansion(mult = c(0, 0.05)),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
  ) +
  labs(x = '', y = '1 - (proportion of miracidia remaining in well)') +
  coord_flip() +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 9),
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 8),
  ) +
  NULL)

save_plot(
  here("Fig5", "plots", "subplots", "Fig5B.pdf"),
  plot,
  base_width = 6,
  base_height = 3
)

######################################################
################### FINAL PLOTTING ###################
######################################################

model <- ggdraw() +
  draw_image(magick::image_read_pdf(here(
    "Fig5",
    "plots",
    "subplots",
    "Fig5A.pdf"
  )))

save_plot(
  here("Fig5", "plots", "Fig5.pdf"),
  plot_grid(
    model,
    plot,
    nrow = 2,
    rel_heights = c(1, 1),
    labels = c('A', 'B'),
    label_y = 1.01
  ),
  base_width = 6,
  base_height = 6
)

save_plot(
  here("Fig5", "plots", "Fig5.png"),
  plot_grid(
    model,
    plot,
    nrow = 2,
    rel_heights = c(1, 1),
    labels = c('A', 'B'),
    label_y = 1.01
  ),
  base_width = 6,
  base_height = 6,
  bg = 'white'
)
