//Stan model to estimate rt seasonal and country variability

// The input data
data {
  int<lower=0> N_countries;
  int<lower=0> N_seasons;
  int<lower=0> N;
  int country[N];
  int season[N];
  vector[N] Rnull;
}

// The parameters accepted by the model
parameters {
  real mean_val;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  vector[N_countries] a;
  vector[N_seasons] b;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  a ~ normal(0, sigma_a);
  b ~ normal(0, sigma_b);
   
  for(i in 1:N){
    log(Rnull[i])  ~  normal(mean_val + a[country[i]] + b[season[i]], sigma);
  }
  
}

generated quantities {
  real country_eff;
  real season_eff;
  real Rnull_country_sim;
  real Rnull_season_sim;
  real Rnull_relative_country_sim;
  real Rnull_relative_season_sim;
  
  country_eff = normal_rng(0,sigma_a);
  Rnull_country_sim = exp(mean_val+country_eff);
  Rnull_relative_country_sim = (exp(mean_val+country_eff)/exp(mean_val) - 1);
  
  season_eff = normal_rng(0,sigma_b);
  Rnull_season_sim = exp(mean_val+season_eff);
  Rnull_relative_season_sim = (exp(mean_val+season_eff)/exp(mean_val) - 1);

}

