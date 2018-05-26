source('R/load-deps.R')

joint_post <- get_joint_post(c('B'))

primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup

translocations <- full_join(
  read_csv('data/out/alpine-translocations.csv') %>%
    mutate(pit_tag_id = as.character(pit_tag_id)),
  read_csv('data/out/subalpine-translocations.csv') %>%
    mutate(pit_tag_id = as.character(pit_tag_id))
)

n_intro_counts <- translocations %>%
  count(translocate_date, site) %>%
  mutate(site = factor(paste(site, 'lake'),
                       levels = c('Alpine', 'Subalpine')),
         group = paste(translocate_date, site),
         label = paste0('n=', n),
         date = as.Date(translocate_date))

intro_df <- n_intro_counts %>%
  left_join(primary_periods) %>%
  select(site, date, n, t) %>%
  mutate(t = t-1)




b_df <- lapply(joint_post, FUN = function(x) {
  reshape2::melt(x, varnames = c('iter', 't')) %>%
    as_tibble
}) %>%
  bind_rows(.id = 'site') %>%
  left_join(intro_df) %>%
  # subtract known introductions
  mutate(n = ifelse(is.na(n), 0, n),
         natB = value - n) %>%
  group_by(iter, site) %>%
  arrange(site, iter, t) %>%
  # compute cumulative recruitment
  mutate(cum_recruits = cumsum(natB)) %>%
  group_by(site, t) %>%
  summarize(med = median(cum_recruits),
            lo = quantile(cum_recruits, .05),
            hi = quantile(cum_recruits, .95)) %>%
  mutate(t = t + 1) %>%
  left_join(primary_periods) %>%
  filter(!is.na(med)) %>%
  ungroup %>%
  mutate(site = factor(site, levels = c('Alpine', 'Subalpine')))

b_df %>%
  ggplot(aes(date, med, color = site)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = site), color = NA,
              alpha = .5) +
  geom_point() +
  geom_line() +
  facet_wrap(~ year, nrow = 1, scales = 'free_x',
             strip.position = 'bottom') +
  xlab('') +
  ylab('Cumulative recruitment') +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_fill_pander('') +
  scale_color_pander('') +
  scale_y_log10(breaks = c(1, 10, 100, 800)) +
  theme(legend.position = 'top')
save_fig('recruitment-time-series', width = 7, height = 3)
