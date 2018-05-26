source('R/load-deps.R')

d <- read_csv('data/alpine_captures_2006-2017.csv') %>%
  separate(period_id, into = c('year_index', 'prim_period', 'sec_period')) %>%
  select(-year_index) %>%
  mutate(pit_tag_id = parse_character(pit_tag_id)) %>%
  left_join(survey_temp_data())

# correct some swabs that were diluted 10x
diluted_swab_ids <- c('RKS18148', 'RKS18152', 'RKS19531', 'RKS21401', 'RKS22223', 'RKS22713')
assert_that(all(diluted_swab_ids %in% d$swab_id))
diluted_idx <- d$swab_id %in% diluted_swab_ids
d$bd_load[diluted_idx] <- d$bd_load[diluted_idx] * 10

translocations <- read_csv('data/alpine_introductions_2006-2017.csv') %>%
  mutate(pit_tag_id = parse_character(parse_number(pit_tag_id)))

dana_snow <- dana_snow_data()

# find number of individuals in augmented data
n_aug <- 1600
M <- compute_M(d, translocations, n_aug)



# find dates and indices for all primary periods
prim_date_df <- d %>%
  mutate(prim_period = paste(survey_year, prim_period, sep = "-")) %>%
  group_by(prim_period) %>%
  summarize(date = min(capture_date),
            J = length(unique(capture_date))) %>%
  mutate(source = 'surveys') %>%
  select(-prim_period) %>%
  full_join(tibble(date = unique(translocations$translocate_date),
                   source = 'translocations',
                   J = 0)) %>%
  full_join(tibble(date = min(.$date) - 7,
                   source = 'imaginary_first_period',
                   J = 0)) %>%
  arrange(date) %>%
  mutate(primary_index = order(date))

dana_snow %>%
  ggplot(aes(date, snow_depth_inches, color = as.factor(format(date, '%Y')))) +
  geom_point() +
  geom_rug(aes(x = as.Date(date)), data = prim_date_df, inherit.aes = FALSE)

# use a bernoulli GAM to learn how snow depth at Dana meadows
# predicts whether there is a survey
snow_d <- prim_date_df %>%
  mutate(date = as.Date(date)) %>%
  right_join(dana_snow) %>%
  mutate(is_survey = as.numeric(!is.na(primary_index)),
         doy = yday(date),
         year = year(date)) %>%
  distinct(date, .keep_all = TRUE) %>%
  filter(year > 2005) %>%
  mutate(snow_depth_inches = ifelse(snow_depth_inches > 150 & is_survey,
                             1,
                             snow_depth_inches))

snow_m <- gam(is_survey ~ s(snow_depth_inches) + s(doy),
              data = snow_d,
              family = binomial())


# while we don't have 10 or more surveys per year, draw from the survey model to
# choose a good date for an augmented primary period
survey_target <- 10

snow_d$is_augmented_prim_period <- snow_d$is_survey
ann_counts <- snow_d %>%
  group_by(year) %>%
  summarize(n_survey = sum(is_augmented_prim_period, na.rm = TRUE))
years_needing_more_periods <- ann_counts$year[ann_counts$n_survey < survey_target]

current_seed <- 1
while(!all(subset(ann_counts, year < 2017)$n_survey >= survey_target)) {
  # predict new surveys
  snow_d$predicted_survey <- predict(snow_m,
                                     newdata = snow_d,
                                     type = 'response')
  set.seed(current_seed)
  snow_d$predicted_survey <- rbinom(nrow(snow_d),
                                    size = 1,
                                    prob = snow_d$predicted_survey)

  # for each prediction, determine whether there is a primary period within 7 or fewer days
  for (i in 1:nrow(snow_d)) {
    if (is.na(snow_d$predicted_survey[i]) | !(snow_d$year[i] %in% years_needing_more_periods)) {
      next
    } else if (snow_d$predicted_survey[i] & !snow_d$is_augmented_prim_period[i]) {
      date_diffs <- abs(snow_d$date[i] -  snow_d$date[snow_d$is_augmented_prim_period == 1])
      any_within_five_days <- any(date_diffs <= 4)
      # if not, prediction becomes augmented primary period
      if (!any_within_five_days) {
        snow_d$is_augmented_prim_period[i] <- 1
      }
    }
  }

  ann_counts <- snow_d %>%
    group_by(year) %>%
    summarize(n_survey = sum(is_augmented_prim_period, na.rm = TRUE))

  years_needing_more_periods <- ann_counts$year[ann_counts$n_survey < survey_target]

  print(data.frame(ann_counts))
  current_seed <- current_seed + 1
}

augmented_primary_period_df <- snow_d %>%
  filter(is_augmented_prim_period == 1,
         date <= max(prim_date_df$date),
         date >= min(prim_date_df$date)) %>%
  mutate()

# find out whether we added too many primary periods forany months
years_to_remove <- augmented_primary_period_df %>%
  group_by(year) %>%
  summarize(n = n(),
            n_actual = sum(is_survey),
            n_added  = n - n_actual) %>%
  mutate(n_to_remove = n - survey_target) %>%
  filter(n_to_remove > 0)

if (nrow(years_to_remove) > 0) {
  # now randomly select which augmented primary period to subtract
  set.seed(1)
  periods_to_subtract <- augmented_primary_period_df %>%
    filter(year %in% years_to_remove$year,
           !is_survey == 1) %>%
    left_join(years_to_remove) %>%
    split(.$year) %>%
    map(~ sample_n(.x, size = unique(.x$n_to_remove))) %>%
    bind_rows()
  primary_periods <- anti_join(augmented_primary_period_df, periods_to_subtract)
} else {
  primary_periods <- augmented_primary_period_df
}

snowpack <- read_csv('data/dana_meadows_swe_2006-2017.csv') %>%
  rename(snow_wc = snow_water_equivalent) %>%
  select(-measure_date)

primary_periods <- primary_periods %>%
  mutate(J = ifelse(is.na(J), 0, J),
         source = ifelse(is.na(source), 'augmented', source)) %>%
  select(-predicted_survey, -is_augmented_prim_period) %>%
  mutate(primary_index = 1:n()) %>%
  group_by(year) %>%
  mutate(is_overwinter = date == max(date)) %>%
  left_join(snowpack) %>%
  ungroup

assert_that(primary_periods %>%
              group_by(year) %>%
              summarize(n = n()) %>%
              `[[`('n') %>%
              max <= survey_target)


T <- nrow(primary_periods)

# find n sec periods per primary period
J <- primary_periods$J

# generate a dataset with secondary period dates and primary period indices
sec_date_df <- d %>%
  mutate(prim_period = paste(survey_year, prim_period, sep = "-")) %>%
  group_by(prim_period) %>%
  mutate(date = min(capture_date),
         J = length(unique(capture_date))) %>%
  distinct(capture_date, .keep_all = TRUE) %>%
  ungroup %>%
  select(capture_date, sec_period, air_temp, date, J) %>%
  left_join(primary_periods) %>%
  rename(primary_date = date)

# verify that primary periods without surveys do not appear in sec_date_df
assert_that(
  !(any(filter(primary_periods, J == 0)$date %in% sec_date_df$capture_date))
  )


# generate state observations
d <- d %>%
  left_join(sec_date_df) %>%
  group_by(primary_index, pit_tag_id) %>%
  mutate(first_cap = min(capture_date),
         is_first_cap = capture_date == first_cap,
         first_swab = !is.na(swab_id) & is_first_cap) %>%
  # discard second swabs if they exist
  mutate(swab_id = ifelse(first_swab, swab_id, NA),
         bd_load = ifelse(first_swab, bd_load, NA)) %>%
  ungroup %>%
  mutate(y = case_when(
              .$bd_load == 0 ~ 2,
              .$bd_load > 0 ~ 3,
              is.na(.$bd_load) & !is.na(.$swab_id) ~ 4
             ),
  swabber = gsub(pattern = '[0-9]+', '', x = swab_id)) %>%
  select(-prim_period)


# produce a data frame with augmented individual capture histories
pit_tags_never_seen <- translocations$pit_tag_id[
  !(translocations$pit_tag_id %in% d$pit_tag_id)]

# check to see that there is only one observed state per individual per prim per
class_check <- d %>%
  group_by(pit_tag_id, primary_index) %>%
  filter(!is.na(y)) %>%
  summarize(n_obs_classes = length(unique(y)),
            classes = paste(unique(y), collapse = ", ")) %>%
  arrange(-n_obs_classes)
# maximum number of observed classes must be 1
assert_that(max(class_check$n_obs_classes) == 1)

# broadcast the observed classes in the detection data (necessary because when
# a swab is not collected on second capture, the above logic will create an NA
# value for the observed class y)
d <- d %>%
  group_by(pit_tag_id, primary_index) %>%
  mutate(y_obs = max(y, na.rm = TRUE),
         n_swab = length(unique(na.omit(swab_id))),
         y_obs = ifelse(n_swab > 0, y_obs, 4)) %>% # if a swab not taken, then
                                        # we only know the individual was alive
  ungroup %>%
  mutate(y = y_obs)


assert_that(sum(pit_tags_never_seen %in% d$pit_tag_id) == 0)
assert_that(sum(is.na(d$y_obs)) == 0)
assert_that(!any(is.na(d$primary_index)))

augmented_df <- tibble(pit_tag_id = c(paste0('aug-', 1:n_aug),
                                      pit_tags_never_seen),
                       y = 1)

y_df <- d %>%
  select(pit_tag_id, primary_index, sec_period, y) %>%
  unite(survey_id, primary_index, sec_period, sep = '_') %>%
  complete(pit_tag_id, survey_id, fill = list(y = 1)) %>%
  separate(survey_id, into = c('primary_period', 'secondary_period')) %>%
  mutate(primary_period = parse_number(primary_period),
         secondary_period = parse_number(secondary_period)) %>%
  arrange(pit_tag_id, primary_period, secondary_period)


# augment y_df with new individuals
y_df <- y_df %>%
  unite(survey_identifier, primary_period, secondary_period) %>%
  full_join(augmented_df) %>%
  complete(pit_tag_id, survey_identifier, fill = list(y = 1)) %>%
  filter(!is.na(survey_identifier)) %>%
  separate(survey_identifier, into = c('primary_period', 'secondary_period')) %>%
  mutate(primary_period = parse_number(primary_period),
         secondary_period = parse_number(secondary_period)) %>%
  arrange(pit_tag_id, primary_period, secondary_period)

assert_that(!any(is.na(y_df$y)))

# augment y_df with y = 0 for all primary periods without any surveys
# with y = 0 representing essentially an NA value!
no_survey_df <- tibble(
  primary_period = filter(primary_periods, J == 0)$primary_index,
  y = 0,
  secondary_period = 1)

y_df <- y_df %>%
  full_join(no_survey_df) %>%
  unite(survey_id, primary_period, secondary_period, sep = "_") %>%
  complete(pit_tag_id, survey_id, fill = list(y = 0)) %>%
  separate(survey_id, into = c('primary_period', 'secondary_period')) %>%
  mutate(primary_period = parse_number(primary_period),
         secondary_period = parse_number(secondary_period)) %>%
  arrange(pit_tag_id, primary_period, secondary_period) %>%
  filter(!is.na(pit_tag_id))

# assert that the only surveys where y = 0 correspond to primary periods with
# no surveys
assert_that(
  all(
    y_df %>%
      filter(y == 0) %>%
      distinct(primary_period) %>%
      unlist == sort(filter(primary_periods, J == 0)$primary_index)
    )
  )


Y <- acast(y_df, pit_tag_id ~ primary_period ~ secondary_period,
           fill = 0,
           value.var = "y")

assert_that(dim(Y)[1] == M)
assert_that(dim(Y)[2] == T)
assert_that(dim(Y)[3] == max(y_df$secondary_period))


# Process data for introductions ------------------------------------------
n_intro <- length(unique(translocations$pit_tag_id))
assert_that(all(translocations$pit_tag_id %in% dimnames(Y)[[1]]))

introduced <- dimnames(Y)[[1]] %in% translocations$pit_tag_id

translocations <- translocations %>%
  mutate(date = as.Date(translocate_date)) %>%
  left_join(primary_periods)

t_intro <- rep(0, M) # primary period index when the animal was introduced
# 0 acts as an NA value
for (i in 1:M) {
  if (introduced[i]) {
    translocate_df_row <- which(translocations$pit_tag_id == dimnames(Y)[[1]][i])
    t_intro[i] <- translocations$primary_index[translocate_df_row]
  }
}

assert_that(mean(t_intro > 0) == mean(introduced))



# Deal with detection data ------------------------------------------------
Jtot <- sum(J)

# create a vector of indices for each secondary period
jvec <- c()
for (i in seq_along(J)) {
  if (J[i] > 0) {
    jvec <- c(jvec, 1:J[i])
  }
}
assert_that(length(jvec) == Jtot)


survey_number_df <- tibble(primary_period = 1:T) %>%
  mutate(secondary_idx = 1)

p_df <- tibble(primary_period = rep(1:T, J),
               secondary_idx = jvec) %>%
  mutate(j_idx = 1:n())


j_idx <- survey_number_df %>%
  dplyr::select(primary_period, secondary_idx) %>%
  full_join(p_df) %>%
  acast(primary_period ~ secondary_idx,
        fill = 0,
        value.var = 'j_idx')

assert_that(nrow(sec_date_df) == Jtot)

detection_covs <- sec_date_df %>%
  arrange(capture_date) %>%
  mutate(j = 1:n(),
         air_temp_c = c(scale(air_temp)),
         air_temp_c = ifelse(is.na(air_temp_c),
                             mean(air_temp_c, na.rm = TRUE),
                             air_temp_c))
# above mutation replaces one NA value with the mean

X_p <- model.matrix( ~ air_temp_c, data = detection_covs)

assert_that(nrow(X_p) == sum(J))

# verify that the primary periods with no surveys have all zeros in j_idx
assert_that(identical(names(which(rowSums(j_idx) == 0)),
                      which(J == 0) %>% as.character))

pos_swab_df <- d %>%
  filter(!is.na(bd_load), bd_load > 0) %>%
  select(pit_tag_id, primary_index, bd_load)


# we need to estimate an M by T matrix that contains the potential loads of each
# individual in each primary period. That is, the Bd load conditional on being
# infected... This matrix is partly observed.
load_observed <- matrix(0, nrow = M, ncol = T)
observed_log_load <- matrix(-9999, nrow = M, ncol = T)
for (i in 1:M) {
  if (dimnames(Y)[[1]][i] %in% pos_swab_df$pit_tag_id) {
    for (prim_period in 1:T) {
      swab_subset <- filter(pos_swab_df, pit_tag_id == dimnames(Y)[[1]][i] &
                              primary_index == prim_period)
      if (nrow(swab_subset) > 0) {
        if (swab_subset$bd_load != 0) {
          print(paste('filling in log loads for i =', i))
          # we are only tracking positive (nonzero) swab values
          load_observed[i, prim_period] <- 1
          observed_log_load[i, prim_period] <- log(swab_subset$bd_load)
        }
      }
    }
  }
}

primary_periods <- primary_periods %>%
  ungroup %>%
  mutate(t_vec = date - min(date))

# Visualize primary periods
primary_periods %>%
  ggplot(aes(date, primary_index, color = source)) +
  geom_point(size = 2) +
  theme_minimal() +
  xlab('Date') +
  ylab('Primary period') +
  scale_color_gdocs('Source')

snow_wc_vec <- primary_periods %>%
  distinct(year, snow_wc) %>%
  arrange(year) %>%
  `[[`('snow_wc')


# Bundle up data for Stan -------------------------------------------------

stan_d <- list(M = M,
               T = T,
               maxJ = max(J),
               J = J,
               Jtot = Jtot,
               Y = Y,
               introduced = introduced,
               t_intro = t_intro,
               X_p = X_p,
               n_p = ncol(X_p),
               j_idx = j_idx,
               any_surveys = ifelse(J > 0, 1, 0),
               n_swab = nrow(pos_swab_df),
               t_swab = pos_swab_df$primary_index,
               log_load = log(pos_swab_df$bd_load / 10000),
               load_observed = load_observed,
               observed_log_load = observed_log_load - log(10000),
               is_overwinter = as.numeric(primary_periods$is_overwinter),
               prim_idx = rep(1:T, J),
               x = snow_wc_vec,
               x_t = c(scale(primary_periods$snow_wc)),
               K = length(unique(d$survey_year)),
               year = primary_periods$year - min(primary_periods$year) + 1)

# assert that there are no missing values in the data
assert_that(!any(lapply(stan_d, function(x) any(is.na(x))) %>% unlist %>% c))

save.image('alpine.RData')
