// A script that fits to previous season, uses the same inition S,I,R values as well as same beta, prop severe, run the model for this year
data {
  int n_week_fit; // number of observed values, weekly
  int n_week_project; // number of projected values, weekly
  array[n_week_fit] int severe_obs_fit; // observed hospitalisations
  array[n_week_fit]int<lower=0,upper=1> severe_obs_notna; // indicating non-missing data with 1, otherwise 0
  real pop; // population size
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
}
transformed data {
  int n_day_fit = 7 * n_week_fit;
  int n_day_project = 7 * n_week_project;
  
  // epi parameters
  real beta = rate_infectious * Rnull;
}
parameters {
  simplex[3] SIR_ini; // S I R
  real<lower=0, upper=1> prop_severe; // proportion of infections that are severe (aka ILIs?)
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for severe obs fit
}
transformed parameters {
  real<lower=0, upper=1> pop_infect;
  real S_init = SIR_ini[1];
  real I_init = SIR_ini[2];
  real R_init = SIR_ini[3];
  real S_ini_logit = logit(S_init);
  real I_ini_logit = logit(I_init);
  real Rnull_eff = Rnull*(1-SIR_ini[3]);
  // time loop
  array[n_day_fit] real<lower=0, upper=1> S;
  array[n_day_fit] real<lower=0, upper=1> I;
  array[n_day_fit] real<lower=0, upper=1> R;
  array[n_day_fit] real<lower=0, upper=1> delta_severe;
  vector[n_day_fit] severe_mean;
  array[n_week_fit] real<lower=0> severe_mean_weekly;
  real phi;
  
  // full wave [from which we "learn" the relevant paramters for the wave in question]
  S[1] = SIR_ini[1];
  I[1] = SIR_ini[2];
  R[1] = SIR_ini[3];
  pop_infect = 0;
  for (t in 2:n_day_fit) {
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
    severe_mean[t] = delta_severe[t]*pop;
    pop_infect = pop_infect + delta_infective_exposures;
  }
  // fill first position in other vectors
  severe_mean[1] = severe_mean[2];
  delta_severe[1] = delta_severe[2];
  
  // convert daily to weekly
  for (i in 1:n_week_fit) {
    int day_start = (i-1)*7+1; 
    int day_end = day_start+6;
    severe_mean_weekly[i] = sum( severe_mean[day_start:day_end] );
  }
  // Overdispersion
  phi = 1 / reciprocal_phi;
  
}
model {
  // starting wave, through scenarios
  for (t in 1:n_week_fit) {
    if (severe_obs_notna[t]==1) severe_obs_fit[t] ~ neg_binomial_2( severe_mean_weekly[t], phi ) ;
  }
}
generated quantities {
  // declare variables
  array[n_day_project] real<lower=0,upper=1> gen_S;
  array[n_day_project] real<lower=0,upper=1> gen_I;
  array[n_day_project] real<lower=0,upper=1> gen_R;
  array[n_day_project] real<lower=0> gen_delta_severe;
  array[n_day_project] real<lower=0> gen_severe_mean;
  array[n_week_project] int<lower=0> gen_severe_obs_project;
  array[n_week_project] int<lower=0> gen_severe_obs_fit;
  array[n_week_project] real gen_severe_mean_weekly;
  //
  real beta_j; 
  beta_j = beta; // Here uncertainty can be added
  
  for (t in 1:n_day_fit) {
    real delta_E;
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    //
    if (t==1) { // set initial conditions
    gen_S[t] = SIR_ini[1];
    gen_I[t] = SIR_ini[2];
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
    gen_delta_severe[t] = delta_infective_exposures * prop_severe; 
    gen_severe_mean[t] = gen_delta_severe[t]*pop ;
    //
    if (t==2) { // also impute the first position
    gen_severe_mean[1] = gen_severe_mean[2];
    gen_delta_severe[1] = gen_delta_severe[2];
    }
    }
  }
  
  // convert daily to weekly
  for (i in 1:n_week_project) {
    int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
    int day_end = day_start+6;
    gen_severe_mean_weekly[i] = sum( gen_severe_mean[day_start:day_end] );
  }
  
  // observation loop: past wave
  for (t in 1:n_week_fit) {
    gen_severe_obs_fit[t] = neg_binomial_2_rng( severe_mean_weekly[t], 1 );
  }
  // observation loop: future wave
  for (t in 1:n_week_project) {
    gen_severe_obs_project[t] = neg_binomial_2_rng(gen_severe_mean_weekly[t], 1 ) ;
  }
  
}
