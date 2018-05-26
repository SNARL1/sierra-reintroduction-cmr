source('R/load-deps.R')

pars <- c('alpha_lambda',
          'beta_lambda_ow',
          'beta_lambda_wc',
          'beta0_lambda',
          'sigma_lambda',
          'alpha_g',
          'beta_g',
          'sigma_g')

joint_post <- get_joint_post(pars)

recruit_par_df <- joint_post %>%
  bind_rows(.id = 'Site') %>%
  gather(var, value, -Site) %>%
  mutate(response = ifelse(grepl('lambda', x = var),
                           'Pr(entry)',
                           'Pr(infected | entry)'),
         Parameter = case_when(
           .$var == 'alpha_lambda' ~ 'Intercept',
           .$var == 'beta_lambda_ow' ~ 'Overwinter effect',
           .$var == 'beta_lambda_wc' ~ 'SWE effect',
           .$var == 'beta0_lambda' ~ 'First\nprimary\nperiod\nadjustment',
           .$var == 'sigma_lambda' ~ 'Annual\nstandard\ndeviation',
           .$var == 'alpha_g' ~ 'Intercept',
           .$var == 'beta_g' ~ 'Bd load effect',
           .$var == 'sigma_g' ~ 'Annual\nstandard\ndeviation'
         ))


recruit_par_df %>%
  ggplot(aes(x = value, y = Parameter, fill = Site)) +
  geom_density_ridges(scale = .95, alpha = .7,
                      rel_min_height = .001) +
  xlab('Parameter value') +
  scale_fill_pander('') +
  ylab('') +
  facet_wrap(~response, scales = 'free') +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'top')
save_fig('recruitment-effects', width = 7, height = 4)


recruit_par_df %>%
  group_by(Site, var, response, Parameter) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  write_csv('out/recruit-pars.csv')
