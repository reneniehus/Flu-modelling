data {
  int n_week; // number of observations, weekly
  array[n_week] int severe_obs; // observed hospitalisations
  real pop; // population size
}
transformed data {
  int n_day = 7 * n_week;
}
parameters {
  real I_ini_logit;
  real prop_severe_logit;
  real S_ini_logit;
  
  //real delta_severe_two_log ; 
}
transformed parameters {
  real<lower=0, upper=1> I_ini=inv_logit(I_ini_logit);
  real<lower=0, upper=1> prop_severe=inv_logit(prop_severe_logit);
  real<lower=0, upper=1> S_ini=inv_logit(S_ini_logit);
  // epi parameters
  real Rnull = 2.5; // https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  real rate_infectious = 0.1042; // 1/7 days
  real beta = rate_infectious * Rnull;
  
  // time loop
  array[n_day] real<lower=0, upper=1> S;
  array[n_day] real<lower=0, upper=1> I;
  array[n_day] real<lower=0, upper=1> R;
  array[n_day] real delta_severe;
  
  //real<lower=0, upper=1> delta_infective_exposures_two ; 
  array[n_day] real<lower=0> severe_mean;
  array[n_week] real<lower=0> severe_mean_weekly;
  // initial conditions: S
  S[1] = S_ini;
  // initial conditions: I
  I[1] = I_ini ;
  // initial conditions: R
  R[1] = 1 - (S[1] + I[1]);
  for (t in 2:n_day) {
    // some local variables (only used in this loop and then forgotten)
    real delta_E;
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    // end: local variables
    delta_infective_exposures = beta*S[t-1] *I[t-1];
    delta_S = -delta_infective_exposures;
    delta_I = delta_infective_exposures - I[t-1]*rate_infectious; 
    delta_R = I[t-1]*rate_infectious; 
    //
    S[t] = S[t-1] + delta_S;
    I[t] = I[t-1] + delta_I;
    R[t] = R[t-1] + delta_R;
    //
    delta_severe[t] = delta_infective_exposures * prop_severe; 
    severe_mean[t] = delta_severe[t]*pop;
  }
  // fill first position in other vectors
  severe_mean[1] = severe_mean[2];
  delta_severe[1] = delta_severe[2];
  
  // convert daily to weekly
  for (i in 1:n_week) {
    int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
    int day_end = day_start+6;
    severe_mean_weekly[i] = sum( severe_mean[day_start:day_end] );
  }
  
}
model {
  // priors
  prop_severe_logit ~ normal( 0,3 );
  S_ini_logit ~ normal( 1.386294,1 ); //
  I_ini_logit ~ normal( -9.21,2 ); //
  // observation process
  for (t in 1:n_week) {
    severe_obs[t] ~ poisson( severe_mean_weekly[t] ) ;
  }
}
generated quantities {
  // declare variables
  real prop_severe_logit_prior;
  real S_ini_logit_prior;
  real I_ini_logit_prior;
  
  real<lower=0,upper=1> gen_beta;
  real<lower=0,upper=1> gen_S_ini;
  
  real<lower=0,upper=1> gen_S[n_day];
  real<lower=0,upper=1> gen_I[n_day];
  real<lower=0,upper=1> gen_R[n_day];
  real<lower=0> gen_delta_severe[n_day];
  array[3,n_day] real<lower=0> gen_severe_mean;
  array[3,n_week] int<lower=0> gen_severe_obs;
  real gen_severe_mean_weekly[3,n_week];
  
  // prior shapes
  prop_severe_logit_prior = normal_rng( 0,3 );
  S_ini_logit_prior = normal_rng( 1.386294,1 );
  I_ini_logit_prior = normal_rng( -9.21,0.5 );
  
  // simulate counterfactual waves
  for (j in 1:3){
    row_vector[3] beta_mod =  [0, +1,   0];
    row_vector[3] S_ini_mod = [0,  0,  -1];
    gen_beta = inv_logit( beta_mod[j]+logit(beta) );
    gen_S_ini = inv_logit( S_ini_mod[j]+logit(S_ini) ) ;
    for (t in 1:n_day) {
    real delta_E;
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    
    //
    if (t==1) { // set initial conditions
    gen_S[t] = gen_S_ini;
    gen_I[t] = I_ini ;
    gen_R[t] = 1 - (gen_S[t] + gen_I[t]);
    } else { // or update
    delta_infective_exposures = gen_beta*gen_S[t-1] *gen_I[t-1];
    delta_S = -delta_infective_exposures;
    delta_I = delta_infective_exposures - gen_I[t-1]*rate_infectious; 
    delta_R = gen_I[t-1]*rate_infectious; 
    //
    gen_S[t] = gen_S[t-1] + delta_S;
    gen_I[t] = gen_I[t-1] + delta_I;
    gen_R[t] = gen_R[t-1] + delta_R;
    //
    gen_delta_severe[t] = delta_infective_exposures * prop_severe; 
    gen_severe_mean[j,t] = gen_delta_severe[t]*pop;
    if (t==2) { // also impute the first position
    gen_severe_mean[j,1] = gen_severe_mean[j,2];
    gen_delta_severe[1] = gen_delta_severe[2];
    }
    }
  }
  }
  
  
  
  // convert daily to weekly
  for (j in 1:3) {
    for (i in 1:n_week) {
    int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
    int day_end = day_start+6;
    gen_severe_mean_weekly[j,i] = sum( gen_severe_mean[j,day_start:day_end] );
  }
  }
  
  // observation loop
  for (j in 1:3) {
    for (t in 1:n_week) {
    gen_severe_obs[j,t] = neg_binomial_2_rng(gen_severe_mean_weekly[j,t],50) ;
  }
  }
  
}
