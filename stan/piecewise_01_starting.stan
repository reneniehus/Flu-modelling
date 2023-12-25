// A script that fits straight lines through groups of points
data {
  int<lower=1> n; // number of observations, weekly
  int<lower=1> n_predict;
  int<lower=1> n_with_predict;
  int<lower=1,upper=n> n_group; // number of groups for the linear fits
  int<lower=1,upper=n_group> group[n]; // group identity of observation
  
  real group_intercept[n_group];
  
  real x_linear[n_with_predict] ; //
  real y[n]; // obervations
  
  real<lower=0> prior_intercept_sd; // prior
  real<lower=0> prior_slope_diff; // prior
  
  real expectation_mu;
  real expectation_sd;
}
transformed data {
  int n_slope_diff = n_group - 1 ;
}
parameters {
  real intercept[n_group]; //
  real slope[n_group]; //
  
  real<lower=0> sigma; // residual error
  real<lower=0> prior_slope_diff_sd; //
  //
   
  
  // real slope_pred = slope[n_group] + normal_rng(0,prior_slope_diff_sd);
}
transformed parameters {
  real mu[n] ;
  real slope_diff[n_slope_diff] ;
  real gen_y[n_with_predict];
  // piece-wise linear expectations
  for (i in 1:n) {
    mu[i] = intercept[ group[i] ] + // intercept
    slope[ group[i] ] * x_linear[i] ; // slope
  }
  
  // differences between adjencent slopes
  for (i in 2:n_group){
    slope_diff[i-1] = slope[i] - slope[i-1] ;
  }
  //
   for (i in 1:n_with_predict) {
    // where data exists
    if (i<=n) {gen_y[i] = mu[i]; }
    // beyond the data
    if (i>n) {
      real intercept_pred = mu[n];
      gen_y[i] = intercept_pred + 
      slope[n_group] * x_linear[i] ; 
    }
  }
  
}
model {
  
  //
  for (i in 1:n) {
    y[i] ~ normal( mu[i] , sigma ) ;
  }
  
  
  // priors (classic)
  for (i in 1:n_group) intercept[i] ~ normal(group_intercept[i], prior_intercept_sd ) ;
  // for (i in 1:n_group) slope[i] ~ normal(0, prior_slope_sd ) ;
  
  // priors (regularisation)
  for (i in 2:n_group) slope_diff[i-1] ~ normal(0 , prior_slope_diff_sd );
  prior_slope_diff_sd ~ exponential( 1/prior_slope_diff );
  // regularisation based on expectatoin (weighting makes softer regularisation)
  for (i in 1:n_with_predict) target += 0.5 * normal_lpdf(gen_y[i] | expectation_mu,expectation_sd);
}
generated quantities {
  // declare generated variables ("gen_variables")
  real gen_y_obs[n_with_predict];
  
  // calculate "gen_variables"
  for (i in 1:n_with_predict) gen_y_obs[i] = gen_y[i] + normal_rng(0,sigma) ;
}
