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
  matrix[n_day_fit, n_age_groups] delta_vax; // daily fraction of newly vaccinated individuals per age group
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
  real ve_inf; // vaccine efficacy on infectiousness
  real ve_susc; // vaccine efficacy on susceptability
  real ve_severe; // vaccine efficacy on infectiousness
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
  // SIR compartments unvaccinated and vaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_u;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_u;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_u;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_v;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_v;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_v;
  array[n_day_fit,n_age_groups] real<lower=0, upper=1> delta_severe;
  array[n_day_fit,n_age_groups] real<lower=0> severe_mean;
  real phi;
  
  // weekly stuff
  array[n_week_fit,n_age_groups] real<lower=0> severe_mean_weekly;
 
  // loop through all days
  for (t in 1:n_day_fit){
    // some local variables (only used in this loop and then forgotten, cannot be constrained)
    real delta_S_u;
    real delta_I_u;
    real delta_R_u;
    real delta_S_v;
    real delta_I_v;
    real delta_R_v;
    real delta_infective_exposures_u;
    real delta_infective_exposures_v;
    
    // end: local variables
    
    if ( season_start[t]==1 ){
      //
      // initiate the compartments based on current season\
      // S I R initial values age dist corrected
      for(a in 1:n_age_groups){
        S_u[t,a] = SIR_ini[season_id[t], a,1] * pop_age_group[a,1] / pop; // rescaled
        I_u[t,a] = SIR_ini[season_id[t], a,2] * pop_age_group[a,1] / pop;
        R_u[t,a] = SIR_ini[season_id[t], a,3] * pop_age_group[a,1] / pop;
        S_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        I_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        R_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        
      }
      
      
    } else {
      for(a in 1:n_age_groups){  
        delta_infective_exposures_u = beta * S_u[t-1,a]  * sum(contact_matrix[ : , a]' .* (I_u[t-1,] + I_v[t-1,]*(1-ve_inf)));
        delta_infective_exposures_v = beta * S_v[t-1,a]  * sum(contact_matrix[ : , a]' .* (I_u[t-1,] + I_v[t-1,]*(1-ve_inf))) * (1 - ve_susc);

        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - I_u[t-1,a] * rate_infectious; 
        delta_I_v = delta_infective_exposures_v - I_v[t-1,a] * rate_infectious; 
        delta_R_u = I_u[t-1,a]*rate_infectious; 
        delta_R_v = I_v[t-1,a]*rate_infectious; 
        
        //
        S_u[t,a] = S_u[t-1,a] + delta_S_u - delta_vax[t-1,a] * S_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        S_v[t,a] = S_v[t-1,a] + delta_S_v + delta_vax[t-1,a] * S_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        I_u[t,a] = I_u[t-1,a] + delta_I_u;
        I_v[t,a] = I_v[t-1,a] + delta_I_v;
        R_u[t,a] = R_u[t-1,a] + delta_R_u - delta_vax[t-1,a] * R_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        R_v[t,a] = R_v[t-1,a] + delta_R_v + delta_vax[t-1,a] * R_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);

        //
        delta_severe[t,a] = (delta_infective_exposures_u + delta_infective_exposures_v * (1-ve_severe)) * prop_severe[season_id[t], a] ; 
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
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S_v;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I_v;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R_v;
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_severe;
  matrix<lower=0>[n_day_project, n_age_groups] gen_severe_mean;
  array[n_week_project, n_age_groups] int<lower=0> gen_severe_obs_project;
  array[n_week_fit, n_age_groups] int<lower=0> gen_severe_obs_fit;
  array[n_week_project, n_age_groups] real gen_severe_mean_weekly;
  
  
  real Rnull_eff[n_season];
  
  //
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  
  real beta_j;
  beta_j = beta; // Here uncertainty can be added
  
  for (t in 1:n_day_project) {
    for (a in 1:n_age_groups) {
      real delta_S_u;
      real delta_I_u;
      real delta_R_u;
      real delta_S_v;
      real delta_I_v;
      real delta_R_v;
      real delta_infective_exposures_u; // to do: make into a matrix (dim per scenario)
      real delta_infective_exposures_v;
      
      //
      if (t==1) { // set initial conditions
      gen_S_u[t,a] = SIR_ini_mu[a,1];
      gen_I_u[t,a] = SIR_ini_mu[a,2];
      gen_R_u[t,a] = 1 - (gen_S_u[t,a] + gen_I_u[t,a]);
      gen_S_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
      gen_I_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
      gen_R_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
      } else { // or update
      delta_infective_exposures_u = beta_j * gen_S_u[t-1,a] *sum(contact_matrix[ : , a]' .* (gen_I_u[t-1,]+gen_I_v[t-1,]*(1-ve_inf)));
      delta_infective_exposures_v = beta_j * gen_S_v[t-1,a] *sum(contact_matrix[ : , a]' .* (gen_I_u[t-1,]+gen_I_v[t-1,]*(1-ve_inf))) * (1-ve_susc);
      delta_S_u = -delta_infective_exposures_u;
      delta_S_v = -delta_infective_exposures_v;
      delta_I_u = delta_infective_exposures_u - gen_I_u[t-1,a] * rate_infectious; 
      delta_I_v = delta_infective_exposures_v - gen_I_v[t-1,a] * rate_infectious; 
      delta_R_u = gen_I_u[t-1,a]*rate_infectious; 
      delta_R_v = gen_I_v[t-1,a]*rate_infectious; 
      //
      gen_S_u[t,a] = gen_S_u[t-1,a] + delta_S_u - delta_vax[t-1,a] * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
      gen_S_v[t,a] = gen_S_v[t-1,a] + delta_S_v + delta_vax[t-1,a] * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
      gen_I_u[t,a] = gen_I_u[t-1,a] + delta_I_u;
      gen_I_v[t,a] = gen_I_v[t-1,a] + delta_I_v;
      gen_R_u[t,a] = gen_R_u[t-1,a] + delta_R_u - delta_vax[t-1,a] * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]); 
      gen_R_v[t,a] = gen_R_v[t-1,a] + delta_R_v + delta_vax[t-1,a] * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
      //
      gen_delta_severe[t,a] = (delta_infective_exposures_u + delta_infective_exposures_v * (1-ve_severe)) * prop_severe_mu[a]; 
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
