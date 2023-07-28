data {
  int n;
  int n_country;
  int n_model;
  int n_treatment;
  int<lower=1,upper=n_country> country_id[n];
  int<lower=1,upper=n_model> model_id[n];
  int<lower=1,upper=n_treatment> treatment[n];
  vector[n] burden_log; // y outcome
}
parameters {
  vector[n_treatment] g;
  vector[n_treatment] alpha[n_country];
  vector[n_treatment] beta[n_model];
  //
  vector<lower=0>[n_treatment] sigma_country;
  vector<lower=0>[n_treatment] sigma_model;
  corr_matrix[n_treatment] Rho_country;
  corr_matrix[n_treatment] Rho_model;
  //
  real sigma;
}
model {
  // probability model
  vector[n] mu;
  for (i in 1:n) {
    
    mu[i] = g[ treatment[i] ] + 
    alpha[ country_id[i] , treatment[i]] + 
    beta[ model_id[i], treatment[i] ];
    
  }
  // compared to the naive model
  // mu = a_model[model_id] + b_x1[model_id].*factor_x1 + b_x2[model_id].*factor_x2;
  burden_log ~ normal(mu, sigma);
  
  // adaptive priors
  alpha ~ multi_normal(rep_vector(0, n_treatment),
  quad_form_diag(Rho_country, sigma_country));
  beta ~ multi_normal(rep_vector(0, n_treatment),
  quad_form_diag(Rho_model, sigma_model));
  // fixed priors
  g ~ normal(15, 5);
  Rho_country ~ lkj_corr(n_treatment);
  Rho_model ~ lkj_corr(n_treatment);
  sigma_country ~ exponential(1);
  sigma_model ~ exponential(1);
  sigma ~ exponential(1);
}
generated quantities{
  
  real gen_burgen_log[n_treatment];
  for (i in 1:n_treatment){
    gen_burgen_log[i] = normal_rng(g[i],sigma);
  }
  
}
