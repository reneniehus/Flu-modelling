// A model that 
data {
  int<lower=1> n; 
  int<lower=1> n_location;
  int<lower=1,upper=n_location>location_id[n];
  vector[n] season_t; 
  real y[n];
}
transformed data {
  vector[1] season_shift;
  season_shift[1] = 0.25; // fixed shift
  vector[n] just_sin =  sin( 2.0*pi()* (season_t/365.0+season_shift[1] )) ;
}
parameters {
  vector[1] y_base;
  vector[n_location] season_fold_log;
  real season_fold_log_mu;
  real <lower=0> season_fold_log_sd;
  real<lower=0> sigma;
}
transformed parameters {
  vector[n] mu ; 
  vector[n_location] season_amp;
  vector[n_location] season_fold;
  for (i in 1:n_location) season_fold[i] = exp(season_fold_log[i]);
  for (i in 1:n_location) season_amp[i] = (season_fold[i]-1)./(season_fold[i]+1);
  
  
  for (i in 1:n){
    mu[i] = (season_amp[ location_id[i] ]*just_sin[i] + 1)*y_base[1] ; // mu = f( amp,t,location_id )
  }
  
}
model {
  for (i in 1:n){
    target += normal_lpdf( y[i] | mu[i] , sigma ) ; // same as Rt ~ normal(mu, sigma)
  }
  // priors
  sigma ~ exponential(1);
  //y_base ~ normal(1,1);
  // hierarchical 
  season_fold_log ~ normal(season_fold_log_mu , season_fold_log_sd) ;
}
generated quantities {
  // declare generated variables ideally with prefix "gen_"
}
