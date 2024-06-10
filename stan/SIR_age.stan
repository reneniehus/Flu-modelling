// A script that fits to previous season, uses the same inition S,I,R values as well as same beta, prop severe, run the model for this year
data {
  // n_scenarios
  // how do you want delta to change +/- 10%
  int n_week_fit; // number of observed values, weekly
  int n_week_project; // number of projected values, weekly
  int n_age_groups; // number of age
  int severe_obs_fit[n_week_fit, n_age_groups]; // observed hospitalisations
  real pop; // population size
  matrix[n_age_groups,1] pop_age_group; // population size per age group 
  matrix[n_age_groups, n_age_groups] contact_matrix; //contact matrix
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
  simplex[3] SIR_init[n_age_groups]; // S I R not age dist corrected
  vector<lower=0, upper=1>[n_age_groups] prop_severe; // proportion of infections that are severe (aka ILIs?)
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for severe obs fit
}
transformed parameters {
  matrix<lower=0, upper=1>[3,n_age_groups] SIR_ini; 
  // S I R initial values age dist corrected
  for(a in 1:n_age_groups){
    SIR_ini[1,a] = SIR_init[a,1] * pop_age_group[a,1] / pop;
    SIR_ini[2,a] = SIR_init[a,2] * pop_age_group[a,1] / pop;
    SIR_ini[3,a] = SIR_init[a,3] * pop_age_group[a,1] / pop;
  }
  
  // time loop
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R;
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] delta_severe;
  matrix[n_day_fit,n_age_groups] severe_mean;
  matrix<lower=0>[n_week_fit,n_age_groups] severe_mean_weekly;
  real phi;
  
  // full wave [from which we "learn" the relevant paramters for the wave in question]
  for(a in 1:n_age_groups){
    S[1,a] = SIR_ini[1,a];
    I[1,a] = SIR_ini[2,a];
    R[1,a] = SIR_ini[3,a];
  }
  for (t in 2:n_day_fit) {
    for(a in 1:n_age_groups){
    // some local variables (only used in this loop and then forgotten, cannot be constrained)
    real delta_S;
    real delta_I;
    real delta_R;
    real delta_infective_exposures;
    // end: local variables
    delta_infective_exposures = beta * S[t-1,a]  * sum(contact_matrix[ : , a]' .* I[t-1,]);
    delta_S = -delta_infective_exposures;
    delta_I = delta_infective_exposures - I[t-1,a] * rate_infectious; 
    delta_R = I[t-1,a]*rate_infectious; 
    //
    S[t,a] = S[t-1,a] + delta_S;
    I[t,a] = I[t-1,a] + delta_I;
    R[t,a] = R[t-1,a] + delta_R;
    //
    delta_severe[t,a] = delta_infective_exposures * prop_severe[a]; 
    severe_mean[t,a] = delta_severe[t,a] * pop_age_group[a,1];
    }
  }
  // fill first position in other vectors
  severe_mean[1] = severe_mean[2];
  delta_severe[1] = delta_severe[2];
  
  // convert daily to weekly
  for (i in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      int day_start = (i-1)*7+1; 
      int day_end = day_start+6;
      severe_mean_weekly[i,a] = sum( severe_mean[day_start:day_end,a] );
    }
  }
  
  // Overdispersion
  phi = 1 / reciprocal_phi;
  
  
}
model {
  // starting wave, through scenarios
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
    severe_obs_fit[t,a] ~ neg_binomial_2( severe_mean_weekly[t,a], phi ) ;
  }
  }
}
generated quantities {
  // declare variables
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_S;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_I;
  matrix<lower=0,upper=1>[n_day_project, n_age_groups] gen_R;
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_severe;
  matrix<lower=0>[n_day_project, n_age_groups] gen_severe_mean;
  array[n_week_project, n_age_groups] int<lower=0> gen_severe_obs_project;
  array[n_week_project, n_age_groups] int<lower=0> gen_severe_obs_fit;
  array[n_week_project, n_age_groups] real gen_severe_mean_weekly;
  //
  real beta_j; 
  beta_j = beta; // Here uncertainty can be added
  
  real Rnull_eff = Rnull*(1-sum(SIR_ini[3,]));
  
  for (t in 1:n_day_fit) {
    for (a in 1:n_age_groups) {
      real delta_E;
      real delta_S;
      real delta_I;
      real delta_R;
      real delta_infective_exposures; // to do: make into a matrix (dim per scenario)
      //
      if (t==1) { // set initial conditions
        gen_S[t,a] = SIR_ini[1,a];
        gen_I[t,a] = SIR_ini[2,a];
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
        gen_delta_severe[t,a] = delta_infective_exposures * prop_severe[a]; 
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
