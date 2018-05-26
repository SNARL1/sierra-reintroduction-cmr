library(tidyverse)
library(ggthemes)
library(viridis)
library(assertthat)
library(patchwork)

theme_set(theme_minimal())


# Plots to demonstrate primary period augmentation ------------------------
primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup %>%
  mutate(Source = case_when(
           .$source %in% c('augmented', 'imaginary_first_period') ~ 'Augmented',
           .$source == 'surveys' ~ 'Surveys',
           .$source == 'translocations' ~ 'Translocations'
         ),
         is_translocation = source == 'translocations',
         Site = ifelse(site == 'Subalpine',
                       'Subalpine Lake',
                       'Alpine Lake'))

plot_one_site <- function(df) {
  df %>%
    ggplot(aes(doy, t, color = Source, shape = Source, size = Source)) +
    geom_point() +
    xlab('Day of year') +
    ylab('Primary period') +
    scale_size_manual('Origin', values = c(1, 2, 3)) +
    scale_shape_manual('Origin', values = c(1, 19, 19)) +
    scale_color_manual('Origin', values = c('black', 2, 'dodgerblue')) +
    facet_wrap(~ year, scales = 'free_y') +
    scale_x_continuous(breaks = c(150, 200, 250)) +
    scale_y_reverse(breaks = seq(1, 200, by = 2)) +
    theme(legend.position = 'none',
          panel.grid.minor = element_blank())
}


alpine_plot <- primary_periods %>%
  filter(site == 'Alpine') %>%
  plot_one_site +
  ggtitle('A. Primary periods at Alpine Lake')

subalpine_plot <- primary_periods %>%
  filter(site == 'Subalpine') %>%
  plot_one_site +
  ggtitle('B. Primary periods at Subalpine Lake')

p <- alpine_plot / subalpine_plot
ggsave(filename = 'fig/primary-periods.pdf',
       plot = p,
       width = 7, height = 7)
