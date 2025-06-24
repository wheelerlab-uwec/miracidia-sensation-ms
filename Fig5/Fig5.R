library(here)
library(tidyverse)
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
  geom_bar(stat = "summary", fun = "mean") +
  ggbeeswarm::geom_quasirandom(size = 2, shape = 21, color = 'white', fill = 'steelblue') +
  stat_summary(fun.data = "mean_cl_normal", color = 'indianred') +
  ggpubr::stat_compare_means(
    comparisons = comparisons,
    label = 'p.signif',
    method = 't.test',
    tip.length = 0,
    size = 2
  ) +
  scale_x_discrete(
    name = NULL,
    labels = labels
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = '', y = 'Fraction successfully penetrating') +
  coord_flip() +
  theme_cowplot() +
  theme(
    axis.title = element_text(size = 9),
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 8),
  ) +
  NULL)

save_plot(here("Fig5", "plots", "Fig5.pdf"), plot, base_width = 6, base_height = 3)
