data {
  //  comment
  int n; 
  int n_model;
  int model_id[n];
  vector[n] factor_x1;
  vector[n] factor_x2;
  vector[n] burden_log;
}
parameters {
  vector[n_model] a_model;
  vector[n_model] b_x1;
  vector[n_model] b_x2;
  
  real mu_bx1;
  real mu_bx2;
  real mu_a;
  
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_bx1;
  real<lower=0> sigma_bx2;
}
model {
  vector[n] mu;
  
  b_x1 ~ normal(mu_bx1,sigma_bx1);
  b_x2 ~ normal(mu_bx2,sigma_bx2);
  a_model ~ normal( mu_a, sigma_a );
  //
  mu_bx1 ~ normal(0, 5);
  mu_bx2 ~ normal(0, 5);
  mu_a ~ normal(0, 5);
  //
  sigma ~ exponential(1); // FIXME: prior check
  sigma_a ~ exponential(1);
  sigma_bx1 ~ exponential(1);
  sigma_bx2 ~ exponential(1);
  
  mu = a_model[model_id] + b_x1[model_id].*factor_x1 + b_x2[model_id].*factor_x2 ;
  burden_log ~ normal(mu, sigma);
}
