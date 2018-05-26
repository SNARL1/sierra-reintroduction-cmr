source('R/load-deps.R')
dir.create('data/out')


if (!file.exists('alpine.RData')) source('R/clean-alpine-data.R')
load('alpine.RData')
m_init <- stan_model('stan/uncertain-state.stan')

m_fit <- vb(m_init,
            data = stan_d, 
            seed = 12,
            eta = .4, 
            init = 0, 
            adapt_engaged = FALSE)
# check to ensure that the variational approximation did not converge to the
# maximum size of the superpopulation (happens rarely)
hist(rstan::extract(m_fit, pars = 'Nsuper')$Nsuper)
write_rds(m_fit, 'alpine_fit.rds')

# save some key data from alpine
write_csv(primary_periods, 'data/out/alpine-pp.csv')
write_csv(pos_swab_df, 'data/out/alpine-pos-swabs.csv')
write_rds(stan_d, 'data/out/alpine-stan-d.rds')
write_rds(detection_covs %>% mutate(site = 'Alpine'),
          'data/out/alpine-detection-covs.rds')
write_rds(Y, 'data/out/alpine-Y.rds')
write_csv(translocations %>% mutate(site = 'Alpine'),
          'data/out/alpine-translocations.csv')

