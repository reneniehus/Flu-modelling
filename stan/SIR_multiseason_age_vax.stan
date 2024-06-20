// A model that fits to several seasons, is age structured, differentiates vaccine status
data {
  // data relevated for the fit 
  int n_season;      // number of seasons
  int n_week_fit;    // number of observable values, weekly
  int n_day_fit;     // number of obervatble values, daily
  int n_age_groups;  // number of age groups
  int severe_obs_fit[n_week_fit, n_age_groups]; // observed hospitalisations
  array[n_week_fit]int<lower=0,upper=1> severe_obs_notna; // indicating non-missing data with 1, otherwise 0
  array[n_day_fit] int<lower=0,upper=2> season_start; // indicating first week of a season with 1, the second week with 2, otherwise 0
  array[n_day_fit] int<lower=1,upper=n_season> season_id; // indicating which seasn each obervable day belongs to
  real pop; // population size
  matrix[n_age_groups,1] pop_age_group; // population size per age group, requires to be a matrix 
  matrix[n_age_groups, n_age_groups] contact_matrix; //contact matrix
  matrix[n_day_fit, n_age_groups] delta_vax; // daily fraction of newly vaccinated individuals per age group
  // data relevant for projected scenarios
  int n_week_project;// number of projected values, weekly
  int n_day_project;
  int n_scenario;// 
  int axis_transmission[n_scenario];
  int axis_vax[n_scenario];
  matrix[n_day_project, n_age_groups] delta_vax_opti;
  matrix[n_day_project, n_age_groups] delta_vax_pess;
  matrix[n_day_project, n_age_groups] delta_vax_null;
  // epi parameters
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
  real ve_inf; // vaccine effectiveness on infectiousness
  real ve_susc; // vaccine effectiveness on susceptability
  real ve_severe; // vaccine effectiveness on severity, given infection
}

transformed data {
  // epi parameters
  real beta = rate_infectious * Rnull;
}

parameters {
  // note: a simplex of length 3, has 2 free parameters
  // note: simplex[n] X[m,o] creates an m x o sized array of simplex, each of size n
  simplex[3] SIR_ini[n_season, n_age_groups]; // S I R initial values per season
  simplex[3] SIR_ini_mu[n_age_groups]; // overall season mean 
  real<lower=0, upper=1> prop_severe[n_season, n_age_groups]; // proportion of infections that are severe (aka ILIs)
  real<lower=0, upper=1> prop_severe_mu[n_age_groups]; // overall mean over season 
  // dispersion parameters
  real<lower=0> sigma_prop_severe;
  real<lower=0> sigma_s;
  real<lower=0> sigma_i;
  
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for severe obs fit
}

transformed parameters {
  // daily stuff
  // SIR compartments unvaccinated and vaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_u; // unvaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_u;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_u;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_v; // vaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_v;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_v;
  array[n_day_fit,n_age_groups] real<lower=0, upper=1> delta_severe;
  array[n_day_fit,n_age_groups] real<lower=0> severe_mean;
  real phi;
  // weekly stuff
  array[n_week_fit,n_age_groups] real<lower=0> severe_mean_weekly;
  
  // Overdispersion
    phi = 1 / reciprocal_phi; // dispersion parameter: var=mu+reciprocal_phi*mu^2
  
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
        S_u[t,a] = SIR_ini[season_id[t], a, 1] * pop_age_group[a, 1] / pop; // rescaled
        I_u[t,a] = SIR_ini[season_id[t], a, 2] * pop_age_group[a, 1] / pop;
        R_u[t,a] = SIR_ini[season_id[t], a, 3] * pop_age_group[a, 1] / pop;
        S_v[t,a] = 0;  // initaillay noone is vaccinated
        I_v[t,a] = 0; 
        R_v[t,a] = 0; 
        
      }
      
    } else {
      for(a in 1:n_age_groups){  
        delta_infective_exposures_u = beta * S_u[t-1,a]  * sum(contact_matrix[ : , a]' .* ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_inf)) );
        delta_infective_exposures_v = beta * S_v[t-1,a]  * sum(contact_matrix[ : , a]' .* ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_inf)) ) * (1 - ve_susc);
        
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - I_u[t-1,a] * rate_infectious; 
        delta_I_v = delta_infective_exposures_v - I_v[t-1,a] * rate_infectious; 
        delta_R_u = I_u[t-1,a]*rate_infectious; 
        delta_R_v = I_v[t-1,a]*rate_infectious; 
        
        S_u[t,a] = S_u[t-1,a] + delta_S_u - delta_vax[t-1,a] * S_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        S_v[t,a] = S_v[t-1,a] + delta_S_v + delta_vax[t-1,a] * S_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        I_u[t,a] = I_u[t-1,a] + delta_I_u;
        I_v[t,a] = I_v[t-1,a] + delta_I_v;
        R_u[t,a] = R_u[t-1,a] + delta_R_u - delta_vax[t-1,a] * R_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        R_v[t,a] = R_v[t-1,a] + delta_R_v + delta_vax[t-1,a] * R_u[t-1,a] / (S_u[t-1,a] + R_u[t-1,a]);
        
        //
        delta_severe[t,a] = (delta_infective_exposures_u * 1 + delta_infective_exposures_v * (1-ve_severe) ) * prop_severe[season_id[t], a] ; 
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
        // define 2 local variables
        int day_start = (i-1)*7+1; 
        int day_end = day_start+6;
        severe_mean_weekly[i,a] = sum( severe_mean[day_start:day_end,a] );
      }
    }
    
  }
}

model {
  
  // --------------------------------likelihood part
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      if (severe_obs_notna[t]==1) severe_obs_fit[t,a] ~ neg_binomial_2( severe_mean_weekly[t,a], phi ) ;
    }
  }
  
  // --------------------------------prior part
  // prior
  for (a in 1:n_age_groups) {
    logit(prop_severe[,a]) ~ normal( logit(prop_severe_mu[a] ) , sigma_prop_severe );
    logit( SIR_ini[,a,1] ) ~ normal( logit(SIR_ini_mu[a,1] ) , sigma_s );
    logit( SIR_ini[,a,2] ) ~ normal( logit(SIR_ini_mu[a,2] ) , sigma_i );
  }
  // priors on dispersion parameters
  sigma_prop_severe    ~ exponential(5); // parameter is the exponential RATE, so BIG numbers mean LOW mean
  sigma_s              ~ exponential(5); // parameter is the exponential RATE, so BIG numbers mean low mean
  sigma_i              ~ exponential(1); // parameter is the exponential RATE, so BIG numbers mean low mean
  
}

generated quantities {
  // --------------------------------declare generated variables
  
  // we give generated quantities the prefix "gen_"
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R_u;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S_v;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I_v;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R_v;
  // note: stan does not have 3-dimensional matrices, thus opting for arrays or 2-dimensional matrixes
  // note: matrix[n,m] M[o] creates an array of length o, each element contraining an nxm matrix, M CONFUSINGLY has then dimension [o,n,m]
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_severe[n_scenario];
  matrix<lower=0>[n_day_project, n_age_groups] gen_severe_mean[n_scenario];
  array[n_week_fit, n_age_groups] int<lower=0> gen_severe_obs_fit;
  array[n_scenario, n_week_project, n_age_groups ] int<lower=0> gen_severe_obs_project;
  array[n_scenario, n_week_project, n_age_groups] real gen_severe_mean_weekly;
  real Rnull_eff[n_season];
  
  // print("dims(gen_delta_severe): ", dims(gen_delta_severe) ); // 
  // print("dims(gen_severe_mean): ", dims(gen_severe_mean) ); // [6,364,1]
  // print("dims(gen_severe_obs_project): ", dims(gen_severe_obs_project) ); // 
  // print("dims(gen_severe_mean_weekly): ", dims(gen_severe_mean_weekly) ); // 
  
  // --------------------------------simulate some quantities of interest
  // computation of Rnull_eff that is not quite correct due to age-structure
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  // --------------------------------simulate fitted observations
  // simulate fitted observations
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      gen_severe_obs_fit[t,a] = neg_binomial_2_rng( severe_mean_weekly[t,a], phi );
    }
  }
  
  // --------------------------------simulate projected observations
  for (j in 1:n_scenario) {
    // settings for the scenarios
    // define 2 local variables
    real beta_j; // scenario-specific beta
    matrix[n_day_project, n_age_groups] delta_vax_j; // scenario-specific vaccine uptake
    if ( axis_transmission[j]==1 ) beta_j = 0.9*beta; // more var: 0.9*normal_rng( log2( beta ), 0.25 ) // pessimistic transmission
    if ( axis_transmission[j]==2 ) beta_j = 1.1*beta; // more var: 1.1*normal_rng( log2( beta ), 0.25 ) // optimisitc transmission
    if ( axis_vax[j]==1 ) delta_vax_j = delta_vax_opti;
    if ( axis_vax[j]==2 ) delta_vax_j = delta_vax_pess;
    if ( axis_vax[j]==3 ) delta_vax_j = delta_vax_null; 
    
    // main time loop
    for (t in 1:n_day_project) {
      for (a in 1:n_age_groups) {
        // define some local variables
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
        gen_S_u[t,a] = gen_S_u[t-1,a] + delta_S_u - delta_vax_j[t-1,a] * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        gen_S_v[t,a] = gen_S_v[t-1,a] + delta_S_v + delta_vax_j[t-1,a] * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        gen_I_u[t,a] = gen_I_u[t-1,a] + delta_I_u;
        gen_I_v[t,a] = gen_I_v[t-1,a] + delta_I_v;
        gen_R_u[t,a] = gen_R_u[t-1,a] + delta_R_u - delta_vax_j[t-1,a] * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]); 
        gen_R_v[t,a] = gen_R_v[t-1,a] + delta_R_v + delta_vax_j[t-1,a] * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        //
        gen_delta_severe[j,t,a] = (delta_infective_exposures_u + delta_infective_exposures_v * (1-ve_severe)) * prop_severe_mu[a]; 
        gen_severe_mean[ j,t,a] = gen_delta_severe[j,t,a] * pop_age_group[a,1] ;
        //
        if (t==2) { // also impute the first position
        gen_severe_mean[j,1,a] = gen_severe_mean[j,2,a];
        gen_delta_severe[j,1,a] = gen_delta_severe[j,2,a];
        }
        }
      } // n_age_groups loop 
    } // n_day_project loop
  } // n_scenario loop
  
  // convert projections daily into weekly
  for (j in 1:n_scenario) {
    for (i in 1:n_week_project) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
        int day_end = day_start+6;
        gen_severe_mean_weekly[j,i,a] = sum( gen_severe_mean[j,day_start:day_end,a] );
      }
    }
  }
  
  // simulate projections
  for (j in 1:n_scenario) {
    for (t in 1:n_week_project) {
      for (a in 1:n_age_groups) {
        gen_severe_obs_project[j,t,a] = neg_binomial_2_rng( gen_severe_mean_weekly[j,t,a] , phi ) ;
      }
    }
  }
  
}
