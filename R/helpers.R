
# Helper functions to process data ----------------------------------------

dana_snow_data <- function() {
  # reads in Dana meadows snow data
  read_csv('data/dana_meadows_snow_depth_2006-2017.csv')
}

survey_temp_data <- function() {
  read_csv('data/ellery_lake_temperature_2006-2017.csv') %>%
    group_by(date) %>%
    summarize(air_temp = mean(temp_c, na.rm = TRUE)) %>%
    rename(capture_date = date)
}

compute_M <- function(d, translocations, n_aug) {
  # compute total number of indiv. in superpopulation
  n_indiv_obs <- length(unique(d$pit_tag_id))
  n_introduced_but_never_seen <- sum(!(translocations$pit_tag_id %in%
                                         d$pit_tag_id))
  M <- n_indiv_obs + n_aug + n_introduced_but_never_seen
  M
}

post2df <- function(x, varnames) {
  x %>%
    reshape2::melt(varnames = varnames) %>%
    as_tibble
}


get_joint_post <- function(pars_to_extract) {
  alpine_post <- rstan::extract(read_rds('alpine_fit.rds'),
                                pars = pars_to_extract)
  subalpine_post <- rstan::extract(read_rds('subalpine_fit.rds'),
                                   pars = pars_to_extract)
  joint_post <- list(Subalpine = subalpine_post,
                     Alpine = alpine_post)
}


save_fig <- function(name, width, height) {
  path <- file.path('fig', name)
  file_names <- paste(path, 'pdf', sep = '.')
  sapply(file_names, ggsave, width = width, height = height)
}

