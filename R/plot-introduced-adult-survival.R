source('R/load-deps.R')

state_df <- get_joint_post('s') %>%
  lapply(FUN = function(x) {
    x[['s']] %>%
      reshape2::melt(varnames = c('iter', 'indiv_idx', 't')) %>%
      as_tibble
  }) %>%
  bind_rows(.id = 'site')


primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup

joint_y <- lapply(c('data/out/alpine-Y.rds', 'data/out/subalpine-Y.rds'), read_rds)
names(joint_y) <- c('Alpine', 'Subalpine')

translocations <- full_join(
  read_csv('data/out/alpine-translocations.csv') %>%
    mutate(pit_tag_id = as.character(pit_tag_id)),
  read_csv('data/out/subalpine-translocations.csv') %>%
    mutate(pit_tag_id = as.character(pit_tag_id))
)



# create a data frame with pit tag ids and individual ids
pit_tag_df <- tibble(pit_tag_id = c(dimnames(joint_y[['Alpine']])[[1]],
                                    dimnames(joint_y[['Subalpine']])[[1]]),
                     indiv_idx = c(1:dim(joint_y[['Alpine']])[[1]],
                                   1:dim(joint_y[['Subalpine']])[[1]]),
                     site = rep(c('Alpine', 'Subalpine'),
                                times = c(dim(joint_y[['Alpine']])[[1]],
                                          dim(joint_y[['Subalpine']])[[1]]))) %>%
  mutate(is_introduced = pit_tag_id %in% translocations$pit_tag_id)


state_summary <- state_df %>%
  left_join(pit_tag_df) %>%
  filter(is_introduced) %>%
  left_join(translocations) %>%
  group_by(site, t, iter, translocate_date) %>%
  summarize(prop_not_recruited = mean(value == 1),
            prop_uninfected = mean(value == 2),
            prop_infected = mean(value == 3),
            prop_dead = mean(value == 4))

n_intro_counts <- translocations %>%
  count(translocate_date, site) %>%
  mutate(group = paste(translocate_date, site),
         label = paste0('n=', n),
         date = as.Date(translocate_date))

intro_survival <- state_summary %>%
  ungroup %>%
  filter(prop_not_recruited != 1) %>%
  group_by(site, t, translocate_date) %>%
  mutate(prop_alive = prop_uninfected + prop_infected) %>%
  summarize(med = median(prop_alive),
            lo = quantile(prop_alive, .05),
            hi = quantile(prop_alive, .95)) %>%
  ungroup %>%
  left_join(dplyr::select(primary_periods, date, t, doy, is_overwinter, site)) %>%
  filter(!is.na(translocate_date)) %>%
  arrange(site, date, translocate_date) %>%
  mutate(group = paste(translocate_date, site))

intro_survival %>%
  ggplot(aes(x = date, y = med, color = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              color = NA, alpha = .6) +
  geom_line() +
  facet_wrap(~site, ncol = 1) +
  xlab('Date') +
  ylab('Proportion surviving') +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1)) +
  scale_fill_solarized() +
  scale_color_solarized() +
  scale_x_date(date_breaks = '12 months', date_labels = "%b %Y",
               limits = c(as.Date('2006-03-01'), max(primary_periods$date))) +
  geom_text(aes(x = date - 120, y = .97, label = label), data = n_intro_counts,
            inherit.aes = FALSE, size = 3.2)
save_fig('introduced-adult-survival', width = 7, height = 3.5)


intro_survival %>%
  write_csv('out/intro-survival.csv')
