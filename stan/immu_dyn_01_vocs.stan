data {
  int n; // number of observations
  int hosp_obs_log[n]; // observed hospitalisations
  int pop;
}
parameters {
  real beta_logit;
  real prop_hosp_logit;
  real S_ini_logit;
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0, upper=1> beta=inv_logit(beta_logit);
  real<lower=0, upper=1> prop_hosp=inv_logit(prop_hosp_logit);
  real<lower=0, upper=1> S_ini=inv_logit(S_ini_logit);
}
model {
  // time loop
  vector[n] S;
  vector[n] E;
  vector[n] I;
  matrix[n,2] R;
  real delta_E;
  real delta_S;
  real delta_I;
  real delta_R;
  real delta_infective_exposures;
  real delta_to_hospital[n];
  //
  beta_logit ~ normal( -0.847, 0.2 );
  prop_hosp_logit ~ normal(-3.89,0.2);
  S_ini_logit ~ normal(1.386,0.2); //
  sigma ~ exponential(0.02);
  
  // initial condition
  S[1] = S_ini;
  I[1] = (exp(hosp_obs_log[1])/pop) /(0.142*prop_hosp); // fix I_ini based on observed hosps
  E[1] = I[1]*0.142/0.385; // coupling E to I based on rates, just guessing
  R[1,1] = 1 - S[1] - I[1] - E[1];
  
  for (t in 2:n) {
    //
    delta_infective_exposures_1 = beta*S[t-1,1] *I[t-1,1];
    delta_infective_exposures_2 = beta*S[t-1,2] *I[t-1,2];
    //
    delta_S_1 = -delta_infective_exposures_1;
    delta_S_2 = -delta_infective_exposures_2;
    delta_E_1 = delta_infective_exposures_1 - E[t-1,1]*0.385; // 1/2.6 days
    delta_E_2 = delta_infective_exposures_2 - E[t-1,2]*0.385; // 1/2.6 days
    delta_I_1 = E[t-1,1]*0.385 - I[t-1,1]*0.142; // 1/7 days
    delta_I_2 = E[t-1,2]*0.385 - I[t-1,2]*0.142; // 1/7 days
    delta_R_1 = I[t-1,1]*0.142; // 1/7 days
    delta_R_2 = I[t-1,2]*0.142; // 1/7 days
    // new VOC
    if ( t==50 ) {
      // swap over new variant to become old
    }
    //
    S[t] = S[t-1] + delta_S;
    E[t] = E[t-1] + delta_E;
    I[t] = I[t-1] + delta_I;
    R[t,1] = R[t-1,1] + delta_R;
    //
    delta_to_hospital[t] = delta_R * prop_hosp;
  }
  // observation process
  for (t in 2:n) {
    hosp_obs_log[t] ~ normal( log( delta_to_hospital[t]*pop ) ,sigma) ;
  }
  
  
}
