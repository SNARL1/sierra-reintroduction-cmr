source('R/load-deps.R')

joint_post <- get_joint_post(c('alpha_phi_infected',
                               'beta_phi_infected',
                               'beta_phi_inf_ow'))

# # make predictions at both sites for individuals at a certain load
bd_load_seq <- 10^c(0:6)
n_primary <- c(Subalpine = 8, Alpine = 10)
assert_that(all(names(n_primary) == names(joint_post)))
load_survival_df <- list()
for (i in seq_along(n_primary)) {
  preds <- array(dim = c(n_primary[i],
                         length(bd_load_seq),
                         length(joint_post[[i]]$alpha_phi_infected)))
  for (j in seq_along(bd_load_seq)) {
    for (k in 1:n_primary[i]) {
      vals <- joint_post[[i]]$alpha_phi_infected +
        joint_post[[i]]$beta_phi_infected * log(bd_load_seq[j] / 10000) +
        joint_post[[i]]$beta_phi_inf_ow * ifelse(k == n_primary[i], 1, 0)
      preds[k, j, ] <- plogis(vals)
    }
  }
  load_survival_df[[i]] <- reshape2::melt(preds,
                                          varnames = c('primary_period',
                                                       'bd_load_idx',
                                                       'iter')) %>%
    as_tibble %>%
    mutate(site = names(n_primary)[i])
}

load_survival_df %>%
  bind_rows %>%
  group_by(site) %>%
  mutate(max_primary = max(primary_period),
         trans = if_else(primary_period == max_primary,
                             'Overwinter transition',
                             'Within-summer transition'),
         trans = factor(trans,
                        levels = c('Within-summer transition',
                                   'Overwinter transition'))) %>%
  filter(primary_period %in% c(1, max_primary)) %>%
  group_by(bd_load_idx, site, trans) %>%
  summarize(med = median(value),
            lo = quantile(value, .05),
            hi = quantile(value, .95)) %>%
  ungroup %>%
  mutate(bd_load = bd_load_seq[bd_load_idx]) %>%
  arrange(site, trans, bd_load) %>%
  ggplot(aes(bd_load, y = med, color = site, fill = site)) +
  facet_wrap(~trans) +
  geom_line() +
  scale_x_log10() +
  geom_ribbon(aes(ymin = lo, ymax = hi), color = NA, alpha = .5) +
  xlab('Bd load') +
  ylab('Survival probability') +
  scale_color_pander('') +
  scale_fill_pander('') +
  theme(legend.position = 'top',
        panel.grid.minor = element_blank())
save_fig('bd-load-survival', width = 7, height = 3)

