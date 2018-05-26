source('R/load-deps.R')

joint_post <- get_joint_post(c('mu_z'))

primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup

pos_swabs <- read_csv('data/out/alpine-pos-swabs.csv') %>%
  mutate(site = 'Alpine',
         pit_tag_id = parse_character(pit_tag_id)) %>%
  full_join(read_csv('data/out/subalpine-pos-swabs.csv') %>%
              mutate(site = 'Subalpine',
                     pit_tag_id = as.character(pit_tag_id))) %>%
  group_by(site) %>%
  mutate(mean_load = mean(bd_load),
         sd_load = sd(bd_load),
         t = primary_index,
         scaled_load = log(bd_load) - log(10000)) %>%
  left_join(primary_periods)


mu_z_df <- joint_post %>%
  lapply(melt, varnames = c('iter', 'year_idx')) %>%
  bind_rows(.id = 'site') %>%
  as_tibble %>%
  full_join(primary_periods) %>%
  group_by(year, site) %>%
  mutate(value = exp(value + log(10000))) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  left_join(primary_periods)

ggplot(mu_z_df, aes(date, med, color = site, fill = site)) +
  geom_jitter(aes(y = bd_load),
              data = pos_swabs %>%
                ungroup %>%
                mutate(site = factor(site, levels = c('Alpine', 'Subalpine'))),
              height = 0, width = 1, alpha = .2) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              color = NA, alpha = .5) +
  scale_y_log10(breaks = c(10^c(-1:8))) +
  facet_wrap(~ year, strip.position = 'bottom', scales = 'free_x',
             nrow = 1) +
  scale_color_pander('') +
  scale_fill_pander('') +
  ylab('Bd load') +
  theme(legend.position = 'top',
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab('')
save_fig('mean-bd-loads', width = 7, height = 3.5)


# Compute correlation in mean load across sites ---------------------------
post_cor <- joint_post %>%
  lapply(melt, varnames = c('iter', 'year_idx')) %>%
  bind_rows(.id = 'site') %>%
  as_tibble %>%
  full_join(primary_periods) %>%
  distinct(site, iter, value, year) %>%
  spread(site, value) %>%
  na.omit %>%
  split(f = .$iter) %>%
  lapply(FUN = function(df) cor(df$Alpine, df$Subalpine)) %>%
  unlist()

post_cor %>%
  data.frame() %>%
  write_csv('out/bd-load-correlations.csv')



# Get coefficient for winter severity ---------------------------------

beta_z <- get_joint_post('beta_z')

beta_z_df <- beta_z %>%
  lapply(melt, varnames = 'iter') %>%
  bind_rows(.id = 'site') %>%
  as_tibble %>%
  group_by(site) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95))

write_csv(beta_z_df, 'out/beta-z-df.csv')
