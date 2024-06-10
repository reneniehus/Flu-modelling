// A script that fits to previous season, uses the same inition S,I,R values as well as same beta, prop severe, run the model for this year
data {
  // n_scenarios;// eg delta change +/- 10%
  int n_season;
  int n_week_fit; // number of observable values, weekly
  int n_day_fit; // number of obervatble values, daily
  int n_week_project; // number of projected values, weekly
  int n_age_groups; // number of age groups
  int severe_obs_fit[n_week_fit, n_age_groups]; // observed hospitalisations
  array[n_week_fit]int<lower=0,upper=1> severe_obs_notna; // indicating non-missing data with 1, otherwise 0
  array[n_day_fit] int<lower=0,upper=2> season_start; // indicating first week of a season with 1, the second week with 2, otherwise 0
  array[n_day_fit] int<lower=1,upper=n_season> season_id; // indicating which seasn each obervable day belongs to
  real pop; // population size
  matrix[n_age_groups,1] pop_age_group; // population size per age group 
  matrix[n_age_groups, n_age_groups] contact_matrix; //contact matrix
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
}

transformed data {
  int n_day_project = 7 * n_week_project;
  // epi parameters
  real beta = rate_infectious * Rnull;
}

parameters {
  simplex[3] SIR_ini[n_season, n_age_groups]; // S I R
  simplex[3] SIR_ini_mu[n_age_groups];// overall mean over season
  
  real<lower=0, upper=1> prop_severe[n_season, n_age_groups]; // proportion of infections that are severe (aka ILIs)
  real<lower=0, upper=1> prop_severe_mu[n_age_groups]; // overall mean over season 
  
  real<lower=0> sigma_prop_severe;
  real<lower=0> sigma_s;
  real<lower=0> sigma_i;
  
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for severe obs fit
}

transformed parameters {
  
  // daily stuff
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R;
  array[n_day_fit,n_age_groups] real<lower=0, upper=1> delta_severe;
  array[n_day_fit,n_age_groups] real<lower=0> severe_mean;
  real phi;
  
  // weekly stuff
  array[n_week_fit,n_age_groups] real<lower=0> severe_mean_weekly;
  
  // loop through all days
  for (t in 1:n_day_fit){
    // some local variables (only used in this loop and then forgotten, cannot be constrained)
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    
    // end: local variables
    
    if ( season_start[t]==1 ){
      //
      // initiate the compartments based on current season\
      // S I R initial values age dist corrected
      for(a in 1:n_age_groups){
        S[t,a] = SIR_ini[season_id[t], a,1] * pop_age_group[a,1] / pop; // rescaled
        I[t,a] = SIR_ini[season_id[t], a,2] * pop_age_group[a,1] / pop;
        R[t,a] = SIR_ini[season_id[t], a,3] * pop_age_group[a,1] / pop;
      }
      
      
    } else {
      for(a in 1:n_age_groups){  
        delta_infective_exposures = beta * S[t-1,a]  * sum(contact_matrix[ : , a]' .* I[t-1,]);
        delta_S = -delta_infective_exposures;
        delta_I = delta_infective_exposures - I[t-1,a] * rate_infectious; 
        delta_R = I[t-1,a]*rate_infectious; 
        //
        S[t,a] = S[t-1,a] + delta_S;
        I[t,a] = I[t-1,a] + delta_I;
        R[t,a] = R[t-1,a] + delta_R;
        //
        delta_severe[t,a] = delta_infective_exposures * prop_severe[season_id[t], a]; 
        severe_mean[t,a] = delta_severe[t,a] * pop_age_group[a,1];
        
        if (season_start[t]==2){
          // fill first position of the season in other vectors
          severe_mean[t-1,a] = severe_mean[t,a];
          delta_severe[t-1,a] = delta_severe[t,a];
        }
      }
      
      
    } // end of daily loop
    
    // convert daily to weekly
    for (i in 1:n_week_fit) {
      for (a in 1:n_age_groups) {
        int day_start = (i-1)*7+1; 
        int day_end = day_start+6;
        severe_mean_weekly[i,a] = sum( severe_mean[day_start:day_end,a] );
      }
    }
    
    // Overdispersion
    phi = 1 / reciprocal_phi; // dispersion parameter: var=mu+reciprocal_phi*mu^2
  }
}

model {
  // starting wave, through scenarios
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      if (severe_obs_notna[t]==1) severe_obs_fit[t,a] ~ neg_binomial_2( severe_mean_weekly[t,a], phi ) ;
    }
  }
  
  for (a in 1:n_age_groups) {
    logit(prop_severe[,a]) ~ normal( logit(prop_severe_mu[a] ) , sigma_prop_severe );
    logit( SIR_ini[,a,1] ) ~ normal( logit(SIR_ini_mu[a,1] ) , sigma_s );
    logit( SIR_ini[,a,2] ) ~ normal( logit(SIR_ini_mu[a,2] ) , sigma_i );
  }
  
  sigma_prop_severe    ~ exponential(5);
  sigma_s ~ exponential(5);
  sigma_i ~ exponential(1);
}

generated quantities {
  // declare variables
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R;
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_severe;
  matrix<lower=0>[n_day_project, n_age_groups] gen_severe_mean;
  array[n_week_project, n_age_groups] int<lower=0> gen_severe_obs_project;
  array[n_week_fit, n_age_groups] int<lower=0> gen_severe_obs_fit;
  array[n_week_fit, n_age_groups] real gen_severe_mean_weekly;
  
  
  real Rnull_eff[n_season];
  
  //
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  
  real beta_j;
  beta_j = beta; // Here uncertainty can be added
  
  for (t in 1:n_day_project) {
    for (a in 1:n_age_groups) {
      real delta_E;
      real delta_S;
      real delta_I;
      real delta_R;
      real delta_infective_exposures; // to do: make into a matrix (dim per scenario)
      //
      if (t==1) { // set initial conditions
      gen_S[t,a] = SIR_ini_mu[a,1];
      gen_I[t,a] = SIR_ini_mu[a,2];
      gen_R[t,a] = 1 - (gen_S[t,a] + gen_I[t,a]);
      } else { // or update
      delta_infective_exposures = beta_j * gen_S[t-1,a] *sum(contact_matrix[ : , a]' .* gen_I[t-1,]);
      delta_S = -delta_infective_exposures;
      delta_I = delta_infective_exposures - gen_I[t-1,a] * rate_infectious; 
      delta_R = gen_I[t-1,a]*rate_infectious; 
      //
      gen_S[t,a] = gen_S[t-1,a] + delta_S;
      gen_I[t,a] = gen_I[t-1,a] + delta_I;
      gen_R[t,a] = gen_R[t-1,a] + delta_R;
      //
      gen_delta_severe[t,a] = delta_infective_exposures * prop_severe_mu[a]; 
      gen_severe_mean[t,a] = gen_delta_severe[t,a] * pop_age_group[a,1] ;
      //
      if (t==2) { // also impute the first position
      gen_severe_mean[1,a] = gen_severe_mean[2,a];
      gen_delta_severe[1,a] = gen_delta_severe[2,a];
      }
      }
    }
  }
  
  // convert daily to weekly
  for (i in 1:n_week_project) {
    for (a in 1:n_age_groups) {
      int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
      int day_end = day_start+6;
      gen_severe_mean_weekly[i,a] = sum( gen_severe_mean[day_start:day_end,a] );
    }
  }
  
  // observation loop: past wave
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      gen_severe_obs_fit[t,a] = neg_binomial_2_rng( severe_mean_weekly[t,a], phi );
    }
  }
  // observation loop: future wave
  for (t in 1:n_week_project) {
    for (a in 1:n_age_groups) {
      gen_severe_obs_project[t,a] = neg_binomial_2_rng(gen_severe_mean_weekly[t,a], phi ) ;
    }
  }
}
