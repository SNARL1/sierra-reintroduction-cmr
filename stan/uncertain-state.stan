data {
  int<lower = 1> M;                   // augmented sample size
  int<lower = 1> T;                   // # primary periods
  int<lower = 1> maxJ;                // max # of secondary periods
  int<lower = 0, upper = maxJ> J[T];  // # 2ndary periods for each prim. period
  int<lower = 1> Jtot;                // total number of surveys
  int<lower = 1, upper = T> prim_idx[Jtot];

  // observations
  // 0=NA, 1=notseen, 2=alive, uninf, 3=alive, infected, 4=alive, unknown
  int<lower = 0, upper = 4> Y[M, T, maxJ];

  int<lower = 0, upper = 1> introduced[M]; // indicator for whether introduced
  int<lower = 0, upper = T> t_intro[M]; // when individuals introduced

  int<lower = 1> n_p;
  matrix[Jtot, n_p] X_p;

  // index order of surveys (0: NA)
  int<lower = 0, upper = Jtot> j_idx[T, maxJ];
  int<lower = 0, upper = 1> any_surveys[T];

  // load data
  int<lower = 1, upper = M * T> n_swab;
  int<lower = 1, upper = T> t_swab[n_swab];   // primary period swab was collected
  vector[n_swab] log_load;                    // observed loads of infected animals
  int<lower = 0, upper = 1> load_observed[M, T];
  matrix[M, T] observed_log_load;
  vector<lower = 0, upper = 1>[T] is_overwinter;

  int K; // number of years
  int<lower = 1, upper = K> year[T];

  // snow_wc
  vector[K] x;
  vector[T] x_t;
}

transformed data {
  int Tm1;                            // # primary periods - 1
  int<lower = 1, upper = 2> introducedP1[M]; // indicator for whether introduced

  Tm1 = T - 1;
  for (i in 1:M)
    introducedP1[i] = introduced[i] + 1;
}

parameters {
  // recruitment
  real alpha_lambda;
  real beta0_lambda;
  real beta_lambda_ow;
  real beta_lambda_wc;
  real<lower = 0> sigma_lambda;
  vector[K] eps_lambda;

  // uninfected survival
  real alpha_phi_uninfected;
  real beta_phi_ow;
  real beta_phi_wc;
  real<lower = 0> sigma_phi_uninf;
  vector[K] eps_phi_uninf;

  // infected survival
  real alpha_phi_infected;
  real beta_phi_inf_wc;
  real beta_phi_infected;
  real beta_phi_inf_ow;
  real<lower = 0> sigma_phi_inf;
  vector[K] eps_phi_inf;

  // loss of infection
  real alpha_eta1;
  real beta_eta1;
  real<lower = 0> sigma_eta1;
  vector[K] eps_eta1;

  // gain of infection
  real alpha_eta2;
  real beta_eta2;
  real<lower = 0> sigma_eta2;
  vector[K] eps_eta2;

  // conditional probability of being infected, given recruitment
  real alpha_g;
  real beta_g;
  real<lower = 0> sigma_g;
  vector[K] eps_g;

  // detection params
  vector[n_p] beta_p;
  real beta_p_infected;
  real beta_p_bd;

  // bd parameters
  real<lower = 0> sigma_bd_swab;
  real alpha_z;
  real<lower = 0> sigma_z;
  vector[K] mu_z;
  matrix[M, T] eps_load;
  real beta_z;

  real<lower = 0, upper = 1> delta; // p(swab succeeds)
}

transformed parameters {
  vector[Jtot] logit_p;
  vector<lower = 0, upper = 1>[Tm1] lambda;
  vector<lower = 0, upper = 1>[Tm1] eta1;
  vector<lower = 0, upper = 1>[Tm1] eta2;
  vector<lower = 0, upper = 1>[Tm1] g;
  vector<lower = 0, upper = 1>[Tm1] gamma[2];   // entry prob
  matrix<lower = -999> [M, T] z; // lower bound protects against indexing errors
  matrix<lower = 0, upper = 1>[M, T] phi_infected;
  matrix<lower = 0, upper = 1>[M, T] phi_uninfected;
  vector[M] log_lik;

  for (i in 1:M) {
    for (t in 1:T) {
      if (load_observed[i, t]) {
        // missing values are encoded as -9999 and the lower bound on
        // potential_log_load ensures that if a missing value is plugged in
        // Stan will raise an error
        z[i, t] = observed_log_load[i, t];
      } else {
        z[i, t] = mu_z[year[t]] + eps_load[i, t];
      }
      phi_uninfected[i, t] = inv_logit(alpha_phi_uninfected
                                       + beta_phi_ow * is_overwinter[t]
                                       + beta_phi_wc * x_t[t]
                                       + eps_phi_uninf[year[t]]);

      phi_infected[i, t] = inv_logit(
        alpha_phi_infected
            + beta_phi_infected * z[i, t]
            + beta_phi_inf_ow * is_overwinter[t]
            + beta_phi_inf_wc * x_t[t]
            + eps_phi_inf[year[t]]);
    }
  }

  // probability of entering population
  lambda = alpha_lambda
            + beta_lambda_ow * is_overwinter[1:Tm1]
            + beta_lambda_wc * x_t[1:Tm1]
            + eps_lambda[year[1:Tm1]];
  lambda[1] = lambda[1] + beta0_lambda;
  lambda = inv_logit(lambda);

  // probability of recruiting as infected, given entry
  g = inv_logit(alpha_g
                + beta_g * mu_z[year[1:Tm1]]
                + eps_g[year[1:Tm1]]);
  gamma[1] = lambda .* (1 - g);
  gamma[2] = lambda .* g;

  eta1 = inv_logit(alpha_eta1
                   + beta_eta1 * mu_z[year[1:Tm1]]
                   + eps_eta1[year[1:Tm1]]);
  eta2 = inv_logit(alpha_eta2
                   + beta_eta2 * mu_z[year[1:Tm1]]
                   + eps_eta2[year[1:Tm1]]);

  logit_p = X_p * beta_p;


  // generate log likelihoods for observation histories
  {
    real acc[4];
    vector[4] gam[T];
    real ps[4, M, Tm1, 4];
    real po[4, M, Jtot, 4];
    real p_uninf;
    real p_inf;

    // define probs of state S(t+1) | S(t)
    // first index: S(t)
    // second index: individual
    // third index: t
    // fourth index: S(t + 1)
    for (i in 1:M) {
      // fill in shared values
      for (t in 1:Tm1) {
        ps[1, i, t, 4] = 0;       // can't die before being alive
        ps[2, i, t, 1] = 0;       // can't unenter population
        ps[3, i, t, 1] = 0;
        ps[2, i, t, 2] = phi_uninfected[i, t] * (1 - eta2[t]);     // survive/didn't get infected
        ps[2, i, t, 3] = phi_uninfected[i, t] * eta2[t]; // survived/became infected
        ps[2, i, t, 4] = 1 - phi_uninfected[i, t]; // death
        ps[3, i, t, 2] = phi_infected[i, t] * eta1[t]; // survived, lost infection
        ps[3, i, t, 3] = phi_infected[i, t] * (1 - eta1[t]); // survived, didn't lose infection
        ps[3, i, t, 4] = 1 - phi_infected[i, t]; // death
        ps[4, i, t, 1] = 0;       // can't unenter after death
        ps[4, i, t, 2] = 0;       // can't un-die
        ps[4, i, t, 3] = 0;
        ps[4, i, t, 4] = 1;       // once dead, always dead
      }

      if (introduced[i]) {
        // individual has been introduced
        // zero probability of recruiting through t_intro - 2
        for (t in 1:(t_intro[i] - 2)) {
          ps[1, i, t, 1] = 1;
          ps[1, i, t, 2] = 0;
          ps[1, i, t, 3] = 0;
        }

        // timestep before introduction has Pr(recruiting) = 1
        // NOTE: this assumes introduced animals are infected upon introduction
        ps[1, i, t_intro[i] - 1, 1] = 0;
        ps[1, i, t_intro[i] - 1, 2] = 0;
        ps[1, i, t_intro[i] - 1, 3] = 1;

        // to avoid NA values, fill in remaining recruitment probs (though they
        // are irrelevant for the likelihood)
        for (t in t_intro[i]:Tm1) {
          ps[1, i, t, 1] = 1 - gamma[1, t] - gamma[2, t];
          ps[1, i, t, 2] = gamma[1, t];
          ps[1, i, t, 3] = gamma[2, t];
        }

      } else {
        for (t in 1:Tm1) {
          ps[1, i, t, 1] = 1 - gamma[1, t] - gamma[2, t];
          ps[1, i, t, 2] = gamma[1, t];
          ps[1, i, t, 3] = gamma[2, t];
        }
      }
    }

    // observation probabilities
    for (i in 1:M) {
      for (j in 1:Jtot) {
        p_uninf = inv_logit(logit_p[j]);
        p_inf = inv_logit(logit_p[j]
                          + beta_p_infected
                          + beta_p_bd * z[i, prim_idx[j]]);
        po[1, i, j, 1] = 1;
        po[1, i, j, 2] = 0;
        po[1, i, j, 3] = 0;
        po[1, i, j, 4] = 0;

        po[2, i, j, 1] = 1 - p_uninf;
        po[2, i, j, 2] = p_uninf * delta;
        po[2, i, j, 3] = 0;
        po[2, i, j, 4] = p_uninf * (1 - delta);

        po[3, i, j, 1] = 1 - p_inf;
        po[3, i, j, 2] = 0;
        po[3, i, j, 3] = p_inf * delta;
        po[3, i, j, 4] = p_inf * (1 - delta);

        po[4, i, j, 1] = 1;
        po[4, i, j, 2] = 0;
        po[4, i, j, 3] = 0;
        po[4, i, j, 4] = 0;
      }
    }

    // likelihood
    for (i in 1:M) {
      // all individuals are in state 1 in first primary period
      gam[1, 1] = 1;
      gam[1, 2] = 0;
      gam[1, 3] = 0;
      gam[1, 4] = 0;

      for (t in 2:T) { // primary periods
        for (k in 1:4) { // state
          for (kk in 1:4) { // previous state
            acc[kk] = gam[t - 1, kk] * ps[kk, i, t - 1, k];
            if (any_surveys[t]) {
              // only increment the log probability with the likelihood
              // if surveys happened
              // (we could equivalently multiply by 1, implying that if there
              //  is no survey, then we cannot make any observation)
              for (j in 1:J[t]) {
                acc[kk] = acc[kk] * po[k, i, j_idx[t, j], Y[i, t, j]];
              }
            }
          }
          gam[t, k] = sum(acc);
        }
      }
      log_lik[i] = log(sum(gam[T]));
    }
  } // end temporary scope
}

model {

  // priors
  alpha_lambda ~ normal(-5, 1.5);
  beta0_lambda ~ normal(0, 1.5);
  beta_lambda_ow ~ normal(0, 1.5);
  beta_lambda_wc ~ normal(0, 1.5);
  sigma_lambda ~ normal(0, 1.5);
  eps_lambda ~ normal(0, sigma_lambda);

  beta_p ~ normal(0, 1.5);
  beta_p_infected ~ normal(0, 1.5);
  beta_p_bd ~ normal(0, 1.5);

  alpha_phi_uninfected ~ normal(0, 1.5);
  beta_phi_ow ~ normal(0, 1.5);
  sigma_phi_uninf ~ normal(0, 1.5);
  eps_phi_uninf ~ normal(0, sigma_phi_uninf);

  alpha_phi_infected ~ normal(0, 1.5);
  beta_phi_infected ~ normal(0, 1.5);
  beta_phi_inf_ow ~ normal(0, 1.5);
  sigma_phi_inf ~ normal(0, 1.5);
  eps_phi_inf ~ normal(0, sigma_phi_inf);

  alpha_eta1 ~ normal(0, 1.5);
  beta_eta1 ~ normal(0, 1.5);
  sigma_eta1 ~ normal(0, 1.5);
  eps_eta1 ~ normal(0, sigma_eta1);

  alpha_eta2 ~ normal(0, 1.5);
  beta_eta2 ~ normal(0, 1.5);
  sigma_eta2 ~ normal(0, 1.5);
  eps_eta2 ~ normal(0, sigma_eta2);

  alpha_g ~ normal(0, 1.5);
  beta_g ~ normal(0, 1.5);
  sigma_g ~ normal(0, 1.5);
  eps_g ~ normal(0, sigma_g);

  sigma_bd_swab ~ normal(0, 1.5);
  to_vector(eps_load) ~ normal(0, sigma_bd_swab);

  alpha_z ~ normal(0, 1.5);
  beta_z ~ normal(0, 1.5);
  sigma_z ~ normal(0, 1.5);
  mu_z ~ normal(alpha_z + beta_z * x, sigma_z);

  delta ~ beta(9, 1);

  beta_phi_wc ~ normal(0, 1.5);
  beta_phi_inf_wc ~ normal(0, 1.5);

  log_load ~ normal(mu_z[year[t_swab]], sigma_bd_swab);
  target += sum(log_lik);
}

generated quantities {
  int<lower = 1, upper = 4> s[M, T];  // latent state
  int<lower=0> Nsuper;                // Superpopulation size
  int<lower = 0> N_uninf[Tm1];
  int<lower = 0> N_inf[Tm1];
  int<lower=0> N[Tm1];                // Actual population size
  int<lower=0> B[Tm1];                // Number of entries

  {
    real ps[4, M, Tm1, 4];
    for (i in 1:M) {
      // fill in shared values
      for (t in 1:(Tm1)) {
        ps[1, i, t, 4] = 0;       // can't die before being alive
        ps[2, i, t, 1] = 0;       // can't unenter population
        ps[3, i, t, 1] = 0;
        ps[2, i, t, 2] = phi_uninfected[i, t] * (1 - eta2[t]);     // survive/didn't get infected
        ps[2, i, t, 3] = phi_uninfected[i, t] * eta2[t]; // survived/became infected
        ps[2, i, t, 4] = 1 - phi_uninfected[i, t]; // death
        ps[3, i, t, 2] = phi_infected[i, t] * eta1[t]; // survived, lost infection
        ps[3, i, t, 3] = phi_infected[i, t] * (1 - eta1[t]); // survived, didn't lose infection
        ps[3, i, t, 4] = 1 - phi_infected[i, t]; // death
        ps[4, i, t, 1] = 0;       // can't unenter after death
        ps[4, i, t, 2] = 0;       // can't un-die
        ps[4, i, t, 3] = 0;
        ps[4, i, t, 4] = 1;       // once dead, always dead
      }

      if (introduced[i]) {
        // individual has been introduced
        // zero probability of recruiting through t_intro - 2
        for (t in 1:(t_intro[i] - 2)) {
          ps[1, i, t, 1] = 1;
          ps[1, i, t, 2] = 0;
          ps[1, i, t, 3] = 0;
        }

        // timestep before introduction has Pr(recruiting) = 1
        // NOTE: this assumes introduced animals are infected upon introduction
        ps[1, i, t_intro[i] - 1, 1] = 0;
        ps[1, i, t_intro[i] - 1, 2] = 0;
        ps[1, i, t_intro[i] - 1, 3] = 1;

        // to avoid NA values, fill in remaining recruitment probs (though they
        // are irrelevant for the likelihood)
        for (t in t_intro[i]:Tm1) {
          ps[1, i, t, 1] = 1 - gamma[1, t] - gamma[2, t];
          ps[1, i, t, 2] = gamma[1, t];
          ps[1, i, t, 3] = gamma[2, t];
        }
      } else {
        for (t in 1:Tm1) {
          ps[1, i, t, 1] = 1 - gamma[1, t] - gamma[2, t];
          ps[1, i, t, 2] = gamma[1, t];
          ps[1, i, t, 3] = gamma[2, t];
        }
      }
    }

    // simulate discrete state values
    for (i in 1:M) {
      s[i, 1] = 1;
      for (t in 2:T) {
        s[i, t] = categorical_rng(to_vector(ps[s[i, t - 1], i, t - 1, ]));
      }
    }
  } // end temporary scope


  {
    int al_inf[M, Tm1];
    int al_uninf[M, Tm1];
    int al[M, Tm1];
    int d[M, Tm1];
    int alive[M];
    int w[M];

    for (i in 1:M) {
      for (t in 2:T) {
        al_uninf[i, t - 1] = (s[i, t] == 2);
        al_inf[i, t - 1] = (s[i, t] == 3);
        al[i, t - 1] = (s[i, t] == 2) + (s[i, t] == 3);
      }
      for (t in 1:Tm1)
        d[i, t] = (s[i, t] == al[i, t]);
      alive[i] = sum(al[i]);
    }

    for (t in 1:Tm1) {
      N_uninf[t] = sum(al_uninf[, t]);
      N_inf[t] = sum(al_inf[, t]);
      N[t] = sum(al[, t]);
      B[t] = sum(d[, t]);
    }
    for (i in 1:M)
      w[i] = 1 - !alive[i];
    Nsuper = sum(w);
  }
}
