source('R/load-deps.R')

eta1_df <- get_joint_post('eta1')
eta2_df <- get_joint_post('eta2')

primary_periods <- read_csv('data/out/alpine-pp.csv') %>%
  mutate(site = 'Alpine') %>%
  full_join(read_csv('data/out/subalpine-pp.csv') %>% mutate(site = 'Subalpine')) %>%
  rename(t = primary_index) %>%
  group_by(site) %>%
  mutate(year_idx = year - min(year) + 1) %>%
  ungroup


# Estimated time series for transition probabilities ----------------------
eta1_df <- lapply(eta1_df, FUN = function(x) {
  x[['eta1']] %>%
    reshape2::melt(varnames = c('iter', 't')) %>%
    as_tibble
}) %>%
  bind_rows(.id = 'site') %>%
  mutate(param = 'eta1')

eta2_df <- lapply(eta2_df, FUN = function(x) {
  x[['eta2']] %>%
    reshape2::melt(varnames = c('iter', 't')) %>%
    as_tibble
}) %>%
  bind_rows(.id = 'site') %>%
  mutate(param = 'eta2')

eta_df <- full_join(eta1_df, eta2_df) %>%
  group_by(site, t, param) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  left_join(primary_periods) %>%
  ungroup %>%
  mutate(Parameter = ifelse(param == 'eta1',
                            'Loss of infection',
                            'Gain of infection'))

transition_plot <- eta_df %>%
  distinct(year, med, lo, hi, Parameter, site) %>%
  ggplot(aes(year, med, color = site)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = site),
              color = NA, alpha = .3) +
  geom_point(size = .5) +
  geom_line() +
  facet_wrap(~ Parameter, scales = 'free_x') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = 'top') +
  scale_color_pander('') +
  scale_fill_pander('') +
  xlab('') +
  ylab('Probability') +
  scale_x_continuous(breaks = seq(2006, 2017, by = 2)) +
  ggtitle('A')


# Visualize transition probabilities as a function of mean bd load --------
joint_post <- get_joint_post(c('mu_z',
                               'alpha_eta1',
                               'beta_eta1',
                               'alpha_eta2',
                               'beta_eta2'))

# Loss of infection -------------------------------------------------------

# Show estimated probabilities across a range of expected Bd values
ng <- 100
eta1_preds <- expand.grid(mu_z = seq(min(joint_post[['Subalpine']]$mu_z),
                                     max(joint_post[['Subalpine']]$mu_z),
                                     length.out = ng),
                          site = c('Subalpine', 'Alpine'),
                          stringsAsFactors = FALSE) %>%
  as_tibble

eta1_preds$med <- NA
eta1_preds$lo <- NA
eta1_preds$hi <- NA

for (i in 1:nrow(eta1_preds)) {
  vals <- joint_post[[eta1_preds$site[i]]]$alpha_eta1 +
    joint_post[[eta1_preds$site[i]]]$beta_eta1 * eta1_preds$mu_z[i]
  vals <- plogis(vals)

  eta1_preds$med[i] <- median(vals)
  eta1_preds$lo[i] <- quantile(vals, .05)
  eta1_preds$hi[i] <- quantile(vals, .95)
}

eta1_probs <- eta1_preds %>%
  mutate(mu_bd = mu_z + log(10000),
         transition = 'Loss of infection')

eta2_preds <- expand.grid(mu_z = seq(min(joint_post[['Subalpine']]$mu_z),
                                     max(joint_post[['Subalpine']]$mu_z),
                                     length.out = ng),
                          site = c('Subalpine', 'Alpine'),
                          stringsAsFactors = FALSE) %>%
  as_tibble

eta2_preds$med <- NA
eta2_preds$lo <- NA
eta2_preds$hi <- NA

for (i in 1:nrow(eta2_preds)) {
  vals <- joint_post[[eta2_preds$site[i]]]$alpha_eta2 +
    joint_post[[eta2_preds$site[i]]]$beta_eta2 * eta2_preds$mu_z[i]
  vals <- plogis(vals)

  eta2_preds$med[i] <- median(vals)
  eta2_preds$lo[i] <- quantile(vals, .025)
  eta2_preds$hi[i] <- quantile(vals, .975)
}

eta2_probs <- eta2_preds %>%
  mutate(mu_bd = mu_z + log(10000),
         transition = 'Gain of infection')


load_transition_plot <- eta1_probs %>%
  full_join(eta2_probs) %>%
  ggplot(aes(exp(mu_bd), med, fill = site, color = site)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              color = NA, alpha = .3) +
  geom_line() +
  facet_wrap(~ transition) +
  xlab('Expected Bd load among infected individuals') +
  scale_x_log10() +
  ylab('Probability') +
  scale_fill_pander('') +
  scale_color_pander('') +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank()) +
  ggtitle('B')

p <- transition_plot / load_transition_plot
ggsave('fig/transition-plot.pdf', plot = p, width = 7, height = 5)
