//Stan model to estimate rt seasonal and country variability

// The input data
data {
  int<lower=0> N_countries;
  int<lower=0> N_seasons;
  int<lower=0> N;
  int country[N];
  int season[N];
  vector[N] Rnull_eff;
  int N_baseline_seasons ; 
  int<lower=1,upper=N_seasons> baseline_seasons[N_baseline_seasons];
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
    log(Rnull_eff[i])  ~  normal(mean_val + a[country[i]] + b[season[i]], sigma);
  }
  
}

generated quantities {
  real country_eff;
  real season_eff;
  real Rnull_eff_country_sim;
  real Rnull_eff_season_sim;
  real Rnull_eff_relative_country_sim;
  real Rnull_eff_relative_season_sim;
  vector[N_baseline_seasons] seasonbaseline;
  real mean_seasonbaseline; 
  real Rnull_eff_relative_seasonbaseline_sim;
  
  // country variability
  country_eff = normal_rng(0,sigma_a);
  Rnull_eff_country_sim = exp(mean_val+country_eff);
  Rnull_eff_relative_country_sim = (exp(mean_val+country_eff)/exp(mean_val) - 1);
  // season variability
  season_eff = normal_rng(0,sigma_b);
  Rnull_eff_season_sim = exp(mean_val+season_eff);
  Rnull_eff_relative_season_sim = (exp(mean_val+season_eff)/exp(mean_val) - 1);
  // season variability relative to 3-season-baseline
  seasonbaseline = mean_val+b[baseline_seasons];
  mean_seasonbaseline = mean(seasonbaseline);
  Rnull_eff_relative_seasonbaseline_sim = (exp(mean_val+season_eff)/exp(mean_seasonbaseline) - 1);

}

