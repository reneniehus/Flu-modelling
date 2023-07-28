// A script that has only 1 scenario (baseline)
data {
  int n_week_full; // number of observations, weekly
  int n_week_start; // number of observations, weekly
  array[n_week_full] int severe_obs_full; // observed hospitalisations
  array[n_week_start] int severe_obs_start; // observed hospitalisations
  real pop_full; // population size
  real pop_start; // population size
}
transformed data {
  int n_day_full = 7 * n_week_full;
  int n_day_start = 7 * n_week_start;
  
  // epi parameters
  real Rnull = 2.0; // https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  real rate_infectious = 0.2777778; // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3057478/#:~:text=We%20estimated%20the%20household%20serial,based%20on%20similar%20study%20designs.
    real beta = rate_infectious * Rnull;
  
  // scenarios
  int n_scenario=1;
  array[n_scenario] real beta_scen;
  array[n_scenario] real susc_red;
  // Behaviour impact on beta
  beta_scen[1] = beta; // baseline
  // beta_scen[2] = beta + 0.10*(1-beta); // X% less social distance/protection
  // Vaccination impact on S
  susc_red[1] = 0; // baseline
}
parameters {
  simplex[3] SIR_ini; // S I R
  real prop_severe_logit;
  // array[n_scenario] real prop_severe_logit_start;
  array[n_scenario] real I_ini_logit_start;
}
transformed parameters {
  real<lower=0, upper=1> prop_severe=inv_logit(prop_severe_logit);
  // array[n_scenario] real<lower=0, upper=1> prop_severe_start=inv_logit(prop_severe_logit); // Changed here from 'prop_severe_logit_start' to 'prop_severe_logit'
  real<lower=0, upper=1> prop_severe_start = prop_severe;
  array[n_scenario] real<lower=0, upper=1> I_ini_start=inv_logit(I_ini_logit_start);
  real<lower=0, upper=1> pop_infect;
  array[n_scenario] real<lower=0, upper=1> pop_infect_start;
  real S_ini_logit = logit(SIR_ini[1]);
  real I_ini_logit = logit(SIR_ini[2]);
  // time loop
  array[n_day_full] real<lower=0, upper=1> S;
  array[n_day_full] real<lower=0, upper=1> I;
  array[n_day_full] real<lower=0, upper=1> R;
  array[n_day_full] real<lower=0, upper=1> delta_severe;
  vector[n_day_full] severe_mean;
  array[n_week_full] real<lower=0> severe_mean_weekly;
  array[n_scenario,n_week_start] real<lower=0> severe_mean_weekly_start;
  
  // full wave [from which we "learn" the relevant paramters (beta) for the wave in question]
  S[1] = SIR_ini[1];
  I[1] = SIR_ini[2];
  R[1] = SIR_ini[3];
  pop_infect = 0;
  for (t in 2:n_day_full) {
    // some local variables (only used in this loop and then forgotten, cannot be constrained)
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    // end: local variables
    delta_infective_exposures = beta * S[t-1] * I[t-1];
    delta_S = -delta_infective_exposures;
    delta_I = delta_infective_exposures - I[t-1] * rate_infectious; 
    delta_R = I[t-1]*rate_infectious; 
    //
      S[t] = S[t-1] + delta_S;
    I[t] = I[t-1] + delta_I;
    R[t] = R[t-1] + delta_R;
    //
      delta_severe[t] = delta_infective_exposures * prop_severe; 
    severe_mean[t] = delta_severe[t]*pop_full;
    pop_infect = pop_infect + delta_infective_exposures;
  }
  // fill first position in other vectors
  severe_mean[1] = severe_mean[2];
  delta_severe[1] = delta_severe[2];
  
  // make fat (fractured-assimilated transmission) wave 
  if (1==0) {
    real mtarget;
    mtarget = sum(severe_mean);
    severe_mean = severe_mean // pile multiple waves on-top of each other
    + 2 * append_row( rep_vector(0,7) , severe_mean[ 1: (n_day_full-1*7)] ) 
    + 1 * append_row( rep_vector(0,2*7) , severe_mean[ 1: (n_day_full-2*7) ] )
    + 0.5 * append_row( rep_vector(0,3*7) , severe_mean[ 1: (n_day_full-3*7) ] ) 
    + 0.25 * append_row( rep_vector(0,4*7) , severe_mean[ 1: (n_day_full-4*7) ] );
    severe_mean = severe_mean/sum(severe_mean);
    severe_mean = severe_mean*mtarget;
  }
  
  // convert daily to weekly
  for (i in 1:n_week_full) {
    int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
    int day_end = day_start+6;
    severe_mean_weekly[i] = sum( severe_mean[day_start:day_end] );
  }
  
  // starting wave [aka the wave in question]: has its own I_ini_start (and NOT prop_severe_start)
  // but same beta 
  // stuff gets overwritten (S,I,R,etc.)
  
  array[n_day_start] real<lower=0, upper=1> S_scen;
  array[n_day_start] real<lower=0, upper=1> I_scen;
  array[n_day_start] real<lower=0, upper=1> R_scen;
  
  for (j in 1:n_scenario){
    // scenario parameters
    real beta_j;
    real S_ini_scen;
    beta_j = beta_scen[j];
    // S_ini_scen = SIR_ini[1]*(1-susc_red[j]);
    S_ini_scen = SIR_ini[1];
    //
    S_scen[1] = S_ini_scen;
    I_scen[1] = I_ini_start[j];
    R_scen[1] = 1 - (S_scen[1] + I_scen[1]);
    
    pop_infect_start[j] = 0;
    for (t in 2:n_day_start) {
      // some local variables (only used in this loop and then forgotten, cannot be constrained)
      // real delta_E;
      real delta_S;
      real delta_I;
      real delta_R;
      real delta_infective_exposures;
      real vaccination;
      // end: local variables
      delta_infective_exposures = beta_j * S_scen[t-1] * I_scen[t-1];
      delta_S = -delta_infective_exposures;
      delta_I = delta_infective_exposures - I_scen[t-1] * rate_infectious; 
      delta_R = I_scen[t-1] * rate_infectious; 
      //
        S_scen[t] = S_scen[t-1] + delta_S;
      I_scen[t] = I_scen[t-1] + delta_I;
      R_scen[t] = R_scen[t-1] + delta_R;
      //
        delta_severe[t] = delta_infective_exposures * prop_severe_start; 
      severe_mean[t] = delta_severe[t] * pop_full;
      pop_infect_start[j] = pop_infect_start[j] + delta_infective_exposures;
      
      // Add vaccination to the population
      if (t == 1e3){
        vaccination = S_scen[t] * susc_red[j];
        S_scen[t] = S_scen[t] - vaccination;
        R_scen[t] = R_scen[t] + vaccination;
      } 
      
    }
    // fill first position in other vectors
    severe_mean[1] = severe_mean[2]; // Absolute numbers
    delta_severe[1] = delta_severe[2]; // Relative numbers: fraction of the population
    
    // make fat (fractured-assimilated transmission) wave 
    if (1==0) {
      real mtarget;
      mtarget = sum(severe_mean);
      severe_mean = severe_mean // pile multiple waves on-top of each other
      + 2 * append_row( rep_vector(0,7) , severe_mean[ 1: (n_day_full-1*7)] ) 
      + 1 * append_row( rep_vector(0,2*7) , severe_mean[ 1: (n_day_full-2*7) ] )
      + 0.5 * append_row( rep_vector(0,3*7) , severe_mean[ 1: (n_day_full-3*7) ] ) 
      + 0.25 * append_row( rep_vector(0,4*7) , severe_mean[ 1: (n_day_full-4*7) ] );
      severe_mean = severe_mean/sum(severe_mean);
      severe_mean = severe_mean*mtarget;
    }
    
    
    // convert daily to weekly
    for (i in 1:n_week_start) {
      int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
      int day_end = day_start+6;
      severe_mean_weekly_start[j,i] = sum( severe_mean[day_start:day_end] );
    }
  }
  
  
}
model {
  // observation process
  // full wave
  for (t in 1:n_week_full) {
    severe_obs_full[t] ~ poisson( severe_mean_weekly[t] ) ;
  }
  // starting wave, through scenarios
  for (j in 1:n_scenario){
    for (t in 1:n_week_start) {
      severe_obs_start[t] ~ poisson( severe_mean_weekly_start[j,t] ) ;
    }
  }
  // priors
  // prop_severe_logit ~ normal( logit(0.02),abs(logit(0.02)-logit(0.20)) );
  prop_severe_logit ~ normal( logit(0.001),abs(logit(0.02)-logit(0.20)) );
  // prop_severe_logit_start ~ normal( logit(0.02),abs(logit(0.02)-logit(0.20)) );
  
  logit(SIR_ini[1]) ~ normal( logit(0.7), abs(logit(0.7)-logit(0.4)) );
  logit(SIR_ini[2]) ~ normal( -10,5 );
  I_ini_logit_start ~ normal( -10,5 ); //
}
generated quantities {
  // declare variables
  real prop_severe_logit_prior;
  real S_ini_logit_prior;
  real I_ini_logit_prior;
  real<lower=0,upper=1> gen_S_ini;
  
  array[n_day_full] real<lower=0,upper=1> gen_S;
  array[n_day_full] real<lower=0,upper=1> gen_I;
  array[n_day_full] real<lower=0,upper=1> gen_R;
  array[n_day_full] real<lower=0> gen_delta_severe;
  
  array[n_scenario,n_day_full] real<lower=0> gen_severe_mean;
  array[n_scenario,n_week_full] int<lower=0> gen_severe_obs;
  array[n_week_full] int<lower=0> gen_severe_obs_full;
  array[n_scenario,n_week_full] real gen_severe_mean_weekly;
  
  // prior shapes
  // prop_severe_logit_prior = normal_rng( logit(0.02),abs(logit(0.02)-logit(0.20)) );
  prop_severe_logit_prior = normal_rng( logit(0.001),abs(logit(0.02)-logit(0.20)) );
  S_ini_logit_prior = normal_rng( logit(0.7), abs(logit(0.7)-logit(0.4)) );
  I_ini_logit_prior = normal_rng( -10,5 );
  
  // simulate counterfactual waves
  for (j in 1:n_scenario){
    // scenario definitions here
    real beta_j;
    beta_j = beta_scen[j];
    // gen_S_ini = SIR_ini[1]*(1-susc_red[j]) ;
    gen_S_ini = SIR_ini[1];
    //
      
      for (t in 1:n_day_full) {
        real delta_E;
        real delta_S;
        real delta_I;
        real delta_R;
        real delta_infective_exposures;
        real vaccination;
        //
          if (t==1) { // set initial conditions
            gen_S[t] = gen_S_ini;
            gen_I[t] = I_ini_start[j] ;
            gen_R[t] = 1 - (gen_S[t] + gen_I[t]);
          } else { // or update
            delta_infective_exposures = beta_j*gen_S[t-1] *gen_I[t-1];
            delta_S = -delta_infective_exposures;
            delta_I = delta_infective_exposures - gen_I[t-1]*rate_infectious; 
            delta_R = gen_I[t-1]*rate_infectious; 
            //
              gen_S[t] = gen_S[t-1] + delta_S;
            gen_I[t] = gen_I[t-1] + delta_I;
            gen_R[t] = gen_R[t-1] + delta_R;
            //
              gen_delta_severe[t] = delta_infective_exposures * prop_severe_start; 
            gen_severe_mean[j,t] = gen_delta_severe[t]*pop_full ;
            if (t==2) { // also impute the first position
              gen_severe_mean[j,1] = gen_severe_mean[j,2];
              gen_delta_severe[1] = gen_delta_severe[2];
            }
            // Add vaccination to the population
            if (t == 1e3){
              vaccination = gen_S[t] * susc_red[j];
              gen_S[t] = gen_S[t] - vaccination;
              gen_R[t] = gen_R[t] + vaccination;
            }
          }
      }
  }
  
  // FIXME: add FAT wave here
  
  // convert daily to weekly
  for (j in 1:n_scenario) {
    for (i in 1:n_week_full) {
      int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
      int day_end = day_start+6;
      gen_severe_mean_weekly[j,i] = sum( gen_severe_mean[j,day_start:day_end] );
    }
  }
  
  // observation loop: Australia
  for (t in 1:n_week_full) {
    gen_severe_obs_full[t] = neg_binomial_2_rng( severe_mean_weekly[t],10000 );
  }
  // observation loop: EU scenarios
  for (j in 1:n_scenario) {
    for (t in 1:n_week_full) {
      gen_severe_obs[j,t] = neg_binomial_2_rng(gen_severe_mean_weekly[j,t],10000) ;
    }
  }
  
}
