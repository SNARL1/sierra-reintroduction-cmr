source('R/load-deps.R')

pars <- c('alpha_phi_uninfected',
          'beta_phi_ow',
          'beta_phi_wc',
          'sigma_phi_uninf')

joint_post <- get_joint_post(pars)

par_df <- joint_post %>%
  bind_rows(.id = 'Site') %>%
  gather(var, value, -Site) %>%
  mutate(Parameter = case_when(
    .$var == 'alpha_phi_uninfected' ~ 'Intercept',
    .$var == 'beta_phi_ow' ~ 'Overwinter effect',
    .$var == 'beta_phi_wc' ~ 'SWE effect',
    .$var == 'sigma_phi_uninf' ~ 'Annual standard deviation'
  )) %>%
  mutate(group = 'Uninfected adults')

infected_pars <- c('alpha_phi_infected',
                   'beta_phi_infected',
                   'beta_phi_inf_ow',
                   'beta_phi_inf_wc',
                   'sigma_phi_inf')

inf_par_df <- get_joint_post(infected_pars) %>%
  bind_rows(.id = 'Site') %>%
  gather(var, value, -Site) %>%
  mutate(Parameter = case_when(
    .$var == 'alpha_phi_infected' ~ 'Intercept',
    .$var == 'beta_phi_infected' ~ 'Bd load effect',
    .$var == 'beta_phi_inf_ow' ~ 'Overwinter effect',
    .$var == 'beta_phi_inf_wc' ~ 'SWE effect',
    .$var == 'sigma_phi_inf' ~ 'Annual standard deviation'
  )) %>%
  mutate(group = 'Bd infected adults')

par_df %>%
  full_join(inf_par_df) %>%
  ggplot(aes(x = value, y = Parameter, fill = Site)) +
  geom_density_ridges(scale = 1, alpha = .7, rel_min_height = .001) +
  xlab('Parameter value') +
  scale_fill_pander('') +
  ylab('') +
  facet_wrap(~group, nrow = 1, scales = 'free_y') +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'top')
save_fig('survival-effects', width = 7, height = 4)

par_df %>%
  full_join(inf_par_df) %>%
  group_by(Site, Parameter, group) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  write_csv('out/survival-params.csv')
