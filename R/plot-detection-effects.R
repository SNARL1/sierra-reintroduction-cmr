source('R/load-deps.R')

joint_post <- get_joint_post(c('beta_p',
                               'beta_p_infected',
                               'beta_p_bd'))

alpine_d <- read_rds('data/out/alpine-stan-d.rds')
subalpine_d <- read_rds('data/out/subalpine-stan-d.rds')

p_post_d <- joint_post %>%
  as.data.frame() %>%
  as_tibble %>%
  gather(param, value) %>%
  mutate(param = gsub('p\\.', replacement = 'p_', x = param)) %>%
  separate(param, into = c('site', 'par'), sep = '\\.') %>%
  mutate(Parameter = case_when(
    .$par == 'beta_p_1' ~ 'Intercept',
    .$par == 'beta_p_2' ~ 'Air temp.',
    .$par == 'beta_p_infected' ~ 'Infected',
    .$par == 'beta_p_bd' ~ 'Bd load'),
    Parameter = factor(Parameter, levels = unique(Parameter)[c(4, 3, 2, 1)]))

p1 <- p_post_d %>%
  ggplot(aes(value, y = Parameter, fill = site)) +
  scale_fill_pander('') +
  geom_density_ridges(alpha = .7, color = alpha(1, .5),
                      bandwidth = .02,
                      rel_min_height = 0.0001) +
  xlab('Value') +
  ylab('Coefficient') +
  theme(panel.grid.minor.x = element_blank()) +
  ggtitle('B') +
  theme(legend.position = 'none')
p1

# Effect of air temp. plot -------------------
logit_p_df <- get_joint_post('logit_p') %>%
  lapply(FUN = function(x) {x[['logit_p']] %>%
      reshape2::melt(varnames = c('iter', 'j'),
                     value.name = 'logit_p') %>%
      as_tibble}) %>%
  bind_rows(.id = 'site')

logit_p_df <- p_post_d %>%
  group_by(site, par) %>%
  mutate(iter = 1:n()) %>%
  select(-Parameter) %>%
  spread(par, value) %>%
  ungroup %>%
  full_join(logit_p_df)

detection_covs <- full_join(read_rds('data/out/alpine-detection-covs.rds') %>%
                              dplyr::select(capture_date, j, air_temp,
                                            is_survey, air_temp_c, site),
                            read_rds('data/out/subalpine-detection-covs.rds') %>%
                              dplyr::select(capture_date, j, air_temp,
                                            is_survey, air_temp_c, site)) %>%
  rename(site = site)

load_vec <- c(alpine_d$log_load, subalpine_d$log_load)

p_summ <- logit_p_df %>%
  full_join(detection_covs) %>%
  mutate(logit_p_infected = logit_p +
           beta_p_infected +
           beta_p_bd * mean(load_vec)) %>%
  dplyr::select(logit_p, logit_p_infected, iter, j, air_temp, site) %>%
  gather(Group, logit_p, -iter, -j, -air_temp, -site) %>%
  mutate(Group = case_when(.$Group == 'logit_p' ~ 'Uninfected adults',
                           .$Group == 'logit_p_infected' ~ 'Bd infected adults'),
         p = plogis(logit_p)) %>%
  group_by(j, air_temp, site, Group) %>%
  summarize(med = median(p),
            lo = quantile(p, .05),
            hi = quantile(p, .95))

p2 <- p_summ %>%
  ungroup %>%
  ggplot(aes(x = air_temp, y = med, group = site,
             fill = site, color = site)) +
  geom_line(alpha = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = site),
              color = NA, alpha = .5) +
  facet_wrap(~Group, nrow = 1) +
  xlab('Survey air temperature (C)') +
  ylab('Probability of detection') +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_color_pander('') +
  scale_fill_pander('') +
  theme(panel.grid.minor = element_blank()) +
  ggtitle('A') +
  theme(legend.position = 'top')
p2

p <- p2 / p1

ggsave('fig/detection-effects.pdf', plot = p, width = 7, height = 5)


# print relevant quantities for manuscript ------------------

# probability of detection
p_summ %>%
  ungroup %>%
  group_by(site) %>%
  mutate(ref_air_temp = 17, 
         dev_air_temp = abs(air_temp - 17), 
         min_dev = min(dev_air_temp, na.rm = TRUE)) %>%
  filter(dev_air_temp == min_dev) %>%
  ungroup %>%
  distinct(air_temp, site, Group, med, lo, hi) %>%
  write_csv(path = 'out/detection-probs.csv')

# coefficient summaries
p_post_d %>%
  group_by(site, Parameter) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  write_csv(path = 'out/detection-effs.csv')
