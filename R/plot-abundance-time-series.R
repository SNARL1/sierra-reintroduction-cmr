source('R/load-deps.R')

joint_post <- get_joint_post(c('N_uninf', 'N_inf'))

primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup

n_df <- lapply(joint_post, FUN = function(x) {
  as.data.frame(x) %>%
    as_tibble() %>%
    gather(var, value) %>%
    separate(var, into = c('var', 't'), sep = '\\.')
}) %>%
  bind_rows(.id = 'site') %>%
  mutate(class = if_else(var == 'N_uninf',
                         'Uninfected',
                         'Infected')) %>%
  group_by(site, class, t) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  mutate(t = parse_number(t)+1) %>%
  full_join(primary_periods) %>%
  ungroup %>%
  arrange(site, t, class)

transloc_summary <- read_csv('data/out/subalpine-translocations.csv') %>%
  mutate(pit_tag_id = as.character(pit_tag_id)) %>%
  full_join(read_csv('data/out/alpine-translocations.csv') %>%
              mutate(pit_tag_id = as.character(pit_tag_id))) %>%
  rename(t = primary_index) %>%
  group_by(t, doy, date, year, site) %>%
  count %>%
  ungroup %>%
  mutate(label_y = ifelse(site == 'Alpine', 320, 80))

dummy_df <- tibble(site = c('Alpine', 'Subalpine'),
                   med = c(100, 100),
                   date = as.Date('2010-06-15'))

abundance_ts <- n_df %>%
  filter(!is.na(class)) %>%
  ggplot(aes(date, med)) +
  geom_blank(data = dummy_df, aes(y = med), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = class, color = class,
                  group = interaction(site, year, class)),
              color = NA, alpha = .4) +
  geom_line(aes(color = class,
                group = interaction(site, year, class))) +
  geom_point(aes(color = class,
                 group = interaction(site, year, class)),
             size = .7) +
  facet_grid(site ~ year, scales = 'free', switch = 'x') +
  xlab('') +
  theme_minimal() +
  ylab('Adult abundance') +
  scale_color_manual('', values = c('red', 'dodgerblue')) +
  scale_fill_manual('', values = c('red', 'dodgerblue')) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'top',
        panel.spacing.x = unit(0, 'lines')) +
  geom_text_repel(aes(x = date, y = label_y, label = paste0('+', n)),
                  inherit.aes = FALSE, color = 'black',
                  data = transloc_summary, size = 3, direction = 'x') +
  geom_vline(aes(xintercept = date), data = transloc_summary,
             linetype = 'dashed', alpha = .3) +
  scale_x_date(labels = date_format("%b"), date_breaks = '1 month') +
  ggtitle('A')
abundance_ts


# Survival time series ----------------------------------------------------
phi_uninf <- get_joint_post('phi_uninfected') %>%
  lapply(FUN = function(x) {
    x %>%
      reshape2::melt(varnames = c('iter', 'individual', 't')) %>%
      as_tibble
  }) %>%
  bind_rows(.id = 'site') %>%
  group_by(t, site) %>%
  summarize(med = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>%
  mutate(group = 'Uninfected')

phi_inf <- get_joint_post('phi_infected') %>%
  lapply(FUN = function(x) {
    x %>%
      reshape2::melt(varnames = c('iter', 'individual', 't')) %>%
      as_tibble
  }) %>%
  bind_rows(.id = 'site') %>%
  group_by(t, site) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  mutate(group = 'Bd infected')


survival_ts_df <- phi_uninf %>%
  full_join(phi_inf) %>%
  left_join(primary_periods) %>%
  mutate(site_class = paste0(site, ': ', group)) %>%
  group_by(site) %>%
  ungroup %>%
  mutate(site = factor(site, levels = c('Alpine', 'Subalpine'))) %>%
  filter(!is.na(site))


survival_ts <- survival_ts_df %>%
  ggplot(aes(date, med, color = group, fill = group,
             group = interaction(year, site, group))) +
  scale_x_date(labels = date_format("%b"), date_breaks = '1 month') +
  facet_grid(site ~ year, scales = 'free_x', switch = 'x') +
  geom_point(size = .7) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .3, color = NA) +
  geom_line(size = .5) +
  scale_color_manual(values = c('red', 'dodgerblue')) +
  scale_fill_manual(values = c('red', 'dodgerblue')) +
  xlab('') +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'none') +
  ylab('Adult survival probability') +
  theme_minimal() +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle('B')
survival_ts

p <- abundance_ts / survival_ts

ggsave(filename = 'fig/abundance-survival-time-series.pdf', plot = p,
       width = 7, height = 5)


# Save summary of Nsuper --------------------------------------------------
N_super <- get_joint_post('Nsuper')

N_summ <- N_super %>%
  bind_rows(.id = 'site') %>%
  group_by(site) %>%
  summarize(med = median(Nsuper),
            lo = quantile(Nsuper, .05),
            hi = quantile(Nsuper, .95))

write_csv(N_summ, 'out/n_summ.csv')
