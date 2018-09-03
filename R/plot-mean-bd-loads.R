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
  mutate(value = exp(value + log(10000))) 

mu_z_summary <- mu_z_df %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  left_join(primary_periods)

load_plot <- mu_z_summary %>%
  ggplot(aes(date, med, color = site, fill = site)) +
  geom_jitter(aes(y = bd_load),
              data = pos_swabs %>%
                ungroup %>%
                mutate(site = factor(site, levels = c('Alpine', 'Subalpine'))),
              height = 0, width = 1, alpha = .2) +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              color = NA, alpha = .4) +
  scale_y_log10(breaks = c(10^c(-1:8))) +
  facet_wrap(~ year, strip.position = 'bottom', scales = 'free_x',
             nrow = 1) +
  scale_color_pander('') +
  scale_fill_pander('') +
  ylab('Bd load') +
  theme(legend.position = 'top',
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(), 
        strip.text = element_blank()) +
  xlab('') + 
  ggtitle('A')
load_plot

# Prevalence plot
prevalence_df <- get_joint_post(c('N_uninf', 'N_inf')) %>%
  lapply(FUN = function(x) {
  as.data.frame(x) %>%
    as_tibble() %>%
    mutate(iter = 1:n()) %>%
    gather(var, value, -iter) %>%
    separate(var, into = c('var', 't'), sep = '\\.')
}) %>%
  bind_rows(.id = 'site') %>%
  mutate(class = if_else(var == 'N_uninf',
                         'Uninfected',
                         'Infected')) %>%
  select(-var) %>%
  spread(class, value) %>%
  group_by(site, iter, t) %>%
  mutate(total = Infected + Uninfected, 
         prevalence = Infected / total) %>%
  filter(total > 0) 

prevalence_summary <- prevalence_df %>%
  group_by(site, t) %>%
  summarize(med = median(prevalence),
            lo = quantile(prevalence, .05),
            hi = quantile(prevalence, .95)) %>%
  mutate(t = parse_number(t)+1) %>%
  full_join(primary_periods) %>%
  ungroup %>%
  arrange(site, t)

prevalence_plot <- prevalence_summary %>%
  ggplot(aes(date, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = site, color = site,
                  group = interaction(site, year)),
              color = NA, alpha = .4) +
  geom_line(aes(color = site,
                group = interaction(site, year))) +
  geom_point(aes(color = site,
                 group = interaction(site, year)),
             size = .7) +
  facet_grid( ~ year, scales = 'free', switch = 'x') +
  xlab('') +
  ylab('Bd prevalence') +
  scale_color_pander('') +
  scale_fill_pander('') +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        panel.spacing.x = unit(0, 'lines')) +
  scale_x_date(labels = date_format("%b"), date_breaks = '1 month') + 
  ggtitle('B')

load_plot / prevalence_plot + plot_layout(heights = c(1, .8))
save_fig('mean-bd-loads', width = 7, height = 4.5)




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




# Compute correlation between mean Bd loads and prevalence ----------------
prevalence_df %>%
  ungroup %>%
  select(site, iter, t, prevalence) %>%
  mutate(t = parse_integer(t)) %>%
  left_join(select(mu_z_df, site, iter, t, value)) %>%
  ungroup %>%
  group_by(site, iter) %>%
  summarize(cor_prev_z = cor(prevalence, value)) %>%
  group_by(site) %>%
  summarize(med = median(cor_prev_z), 
           lo = quantile(cor_prev_z, .05), 
           hi = quantile(cor_prev_z, .95)) %>%
  write_csv('out/cor_prev_z.csv')
