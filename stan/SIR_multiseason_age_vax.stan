// A model that fits to several seasons, is age structured, differentiates vaccine status
data {
  // data relevated for the fit 
  int n_season;      // number of seasons
  int n_week_fit;    // number of observable values, weekly
  int n_day_fit;     // number of obervatble values, daily
  int n_age_groups;  // number of age groups
  int ili_obs_fit[n_week_fit, n_age_groups]; // observed hospitalisations
  array[n_week_fit,n_age_groups]int<lower=0,upper=1> ili_obs_notna; // indicating non-missing data with 1, otherwise 0
  array[n_day_fit] int<lower=0,upper=2> season_start; // indicating first week of a season with 1, the second week with 2, otherwise 0
  array[n_day_fit] int<lower=1,upper=n_season> season_id; // indicating which seasn each obervable day belongs to
  real pop; // population size
  matrix[n_age_groups,1] pop_age_group; // population size per age group, requires to be a matrix 
  matrix[n_age_groups, n_age_groups] contact_matrix; //contact matrix
  matrix[n_day_fit, n_age_groups] delta_vax; // daily fraction of newly vaccinated individuals per age group
  // data relevant for projected scenarios
  int n_week_project; // number of projected weeks
  int n_day_project; // number of projected days
  int n_scenario;// // number of projected scenarios
  int axis_transmission[n_scenario]; // indicator for the transmission scenario axis
  int axis_vax[n_scenario]; // indicator for the vaccine scenario axis
  matrix[n_day_project, n_age_groups] delta_vax_real; // daily assumed vax uptake in projection period
  matrix[n_day_project, n_age_groups] delta_vax_opti; // daily assumed vax uptake in projection period
  matrix[n_day_project, n_age_groups] delta_vax_pess; // daily assumed vax uptake in projection period
  matrix[n_day_project, n_age_groups] delta_vax_null; // daily assumed vax uptake in projection period
  // epi parameters
  real Rnull; // R0
  real rate_infectious; // infectious rate, such that beta = Rnull*rate_infectious
  // (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf)
  real ve_spread; // vaccine effectiveness on onward transmission/infectiousness
  real ve_inf; // vaccine effectiveness on susceptability
  real ve_ili_cond_inf; // vaccine effectiveness on severity, given infection
  //
  int n_daily_time_steps; // number of daily steps
  array[n_day_fit*n_daily_time_steps] int<lower=1,upper=n_day_fit> daily_counter_fit;
}

transformed data {
  // epi parameters
  real beta = rate_infectious * Rnull;
  int n_multi_day_project = n_day_project * n_daily_time_steps;
  int n_multi_day_fit = n_day_fit * n_daily_time_steps;
  real dt = 1.0/n_daily_time_steps;
}

parameters {
  // note: a simplex of length 3, has 2 free parameters
  // note: simplex[n] X[m,o] creates an m x o sized array of simplex, each of size n
  simplex[3] SIR_ini_scaled[n_season, 1]; // S I R initial values per season, 1 can be replaced by n_age_groups
  simplex[3] SIR_ini_mu[1]; // overall season mean 
  real<lower=0, upper=1> prop_ili_scaled[n_season, n_age_groups]; // proportion of infections that are ili 
  real<lower=0, upper=1> prop_ili_mu[n_age_groups]; // overall mean over season 
  // dispersion parameters
  real<lower=0> sigma_prop_ili;
  real<lower=0> sigma_s;
  real<lower=0> sigma_i;
  
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for ili obs fit, var=mu+reciprocal_phi*mu^2
}

transformed parameters {
  // New transformed variables - scaling one element of a simplex or vector
  real<lower=0, upper=1> prop_ili[n_season, n_age_groups];
  array[n_season,1,3] real<lower=0, upper=1> SIR_ini;
  real sum_val;
  for (i in 1:n_season) {
    for (j in 1:n_age_groups) {
      prop_ili[i, j] = prop_ili_scaled[i, j] * 0.1;
    }
    SIR_ini[i,1,1] = SIR_ini_scaled[i,1,1];
    SIR_ini[i,1,2] = SIR_ini_scaled[i,1,2] * 0.01;
    SIR_ini[i,1,3] = SIR_ini_scaled[i,1,3];
    sum_val = SIR_ini[i,1,1] + SIR_ini[i,1,2] + SIR_ini[i,1,3];
    SIR_ini[i,1,1] = SIR_ini[i,1,1] / sum_val;
    SIR_ini[i,1,2] = SIR_ini[i,1,2] / sum_val;
    SIR_ini[i,1,3] = SIR_ini[i,1,3] / sum_val;
  }
  
  // daily stuff
  // SIR compartments unvaccinated and vaccinated
  // variables starting with delta_ are incidence variables that are daily, unless it is specificed otherwise (e.g. through _weekly )
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_u; // susceptible compartment, relative to population size, unvaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_u; // infetious compartment,   relative to population size, unvaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_u; // recovered compartment,   relative to population size, unvaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_v; // susceptible compartment, relative to population size, vaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_v; // infetious compartment,   relative to population size, vaccinated
  matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_v; // recovered compartment,   relative to population size, vaccinated
  array[n_day_fit,n_age_groups] real<lower=0, upper=1> delta_ili; // ili/detectable incidence relative to population size
  array[n_day_fit,n_age_groups] real<lower=0> delta_ili_abs; // ili/detectable incidence in absolute numbers
  real phi; // dispersion parameter, var=mu+reciprocal_phi*mu^2
  // weekly stuff
  array[n_week_fit,n_age_groups] real<lower=0> delta_ili_abs_weekly; // ili/detectable incidence in absolute numbers, weekly aggregate
  
  // Overdispersion
  phi = 1 / reciprocal_phi; // dispersion parameter: var=mu+reciprocal_phi*mu^2
  
  
  
  
  // Declare temporary variables outside the loop
  real prev_S_u[n_age_groups];
  real prev_I_u[n_age_groups];
  real prev_R_u[n_age_groups];
  real prev_S_v[n_age_groups];
  real prev_I_v[n_age_groups];
  real prev_R_v[n_age_groups];
  real prev_delta_ili[n_age_groups];
  real prev_delta_ili_abs[n_age_groups];
  
  real curr_S_u[n_age_groups];
  real curr_I_u[n_age_groups];
  real curr_R_u[n_age_groups];
  real curr_S_v[n_age_groups];
  real curr_I_v[n_age_groups];
  real curr_R_v[n_age_groups];
  real curr_delta_ili[n_age_groups];
  real curr_delta_ili_abs[n_age_groups];
  
  
  // loop through all days
  for (t in 1:n_multi_day_fit){
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
    
    // If it is first day of season AND first iteration of the day
    if ( season_start[daily_counter_fit[t]]==1 && ((t-1)%n_daily_time_steps)==0 ){
      //
      // initiate the compartments based on current season\
      // S I R initial values age dist corrected
      for(a in 1:n_age_groups){
        curr_S_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 1] * pop_age_group[a, 1] / pop; // rescaling
        curr_I_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 2] * pop_age_group[a, 1] / pop; // rescaling
        curr_R_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 3] * pop_age_group[a, 1] / pop; // rescaling
        curr_S_v[a] = 0;  // at start of season, no one is vaccinated
        curr_I_v[a] = 0;
        curr_R_v[a] = 0;
        curr_delta_ili[a] = 1;
        curr_delta_ili_abs[a] = 1;
        
        if (t % n_daily_time_steps == 0) { // Only store every N-th iteration
        S_u[t/n_daily_time_steps,a] = curr_S_u[a];
        I_u[t/n_daily_time_steps,a] = curr_I_u[a];
        R_u[t/n_daily_time_steps,a] = curr_R_u[a];
        S_v[t/n_daily_time_steps,a] = curr_S_v[a];
        I_v[t/n_daily_time_steps,a] = curr_I_v[a];
        R_v[t/n_daily_time_steps,a] = curr_R_v[a];
        delta_ili[t/n_daily_time_steps,a] = curr_delta_ili[a];
        delta_ili_abs[t/n_daily_time_steps,a] = curr_delta_ili_abs[a];
        }
        
      }
      
    } else {
      for(a in 1:n_age_groups){
        //
        delta_infective_exposures_u = dt * beta * prev_S_u[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta * prev_S_v[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        //
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - prev_I_v[a] * rate_infectious*dt;
        delta_R_u = prev_I_u[a]*rate_infectious*dt;
        delta_R_v = prev_I_v[a]*rate_infectious*dt;
        //
        curr_S_u[a] = prev_S_u[a] + delta_S_u - (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_S_v[a] = prev_S_v[a] + delta_S_v + (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_I_u[a] = prev_I_u[a] + delta_I_u;
        curr_I_v[a] = prev_I_v[a] + delta_I_v;
        curr_R_u[a] = prev_R_u[a] + delta_R_u - (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_R_v[a] = prev_R_v[a] + delta_R_v + (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        //
        curr_delta_ili[a] = (delta_infective_exposures_u * 1 + delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili[season_id[daily_counter_fit[t]], a];
        curr_delta_ili_abs[a] = curr_delta_ili[a] * pop_age_group[a,1];
        
        // Only store every n_daily_time_steps-th iteration
        if (t % n_daily_time_steps == 0) { 
        S_u[t/n_daily_time_steps,a] = curr_S_u[a];
        S_v[t/n_daily_time_steps,a] = curr_S_v[a];
        I_u[t/n_daily_time_steps,a] = curr_I_u[a];
        I_v[t/n_daily_time_steps,a] = curr_I_v[a];
        R_u[t/n_daily_time_steps,a] = curr_R_u[a];
        R_v[t/n_daily_time_steps,a] = curr_R_v[a];
        delta_ili[t/n_daily_time_steps,a] = curr_delta_ili[a];
        delta_ili_abs[t/n_daily_time_steps,a] = curr_delta_ili_abs[a];
        }
      }
    }
    // Update previous values with current values for next iteration
    prev_S_u = curr_S_u;
    prev_I_u = curr_I_u;
    prev_R_u = curr_R_u;
    prev_S_v = curr_S_v;
    prev_I_v = curr_I_v;
    prev_R_v = curr_R_v;
    prev_delta_ili = curr_delta_ili;
    prev_delta_ili_abs = curr_delta_ili_abs;
  } // end of multi-daily loop
  
  // convert daily to weekly
  for (i in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      // define 2 local variables
      int day_start = (i-1)*7+1;
      int day_end = day_start+6;
      delta_ili_abs_weekly[i,a] = sum( delta_ili_abs[day_start:day_end,a] );
    }
  }
}

model {
  
  // --------------------------------likelihood part
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      if (ili_obs_notna[t,a]==1) ili_obs_fit[t,a] ~ neg_binomial_2( delta_ili_abs_weekly[t,a]+1e-6, phi ) ; // TODO: remove this line and see if you get priors back
    }
  }
  
  // --------------------------------prior part
  // prior
  for (a in 1:n_age_groups) {
    prop_ili_scaled[,a] ~ uniform(0, 1);
  }
  
  // Uniform priors for SIR_ini
  SIR_ini_scaled[,1,2]   ~ uniform( 0, 1 );
  
  // more priors: put priors on things that we want to fixate more
  logit(prop_ili_mu) ~ normal( logit(0.015) , 0.5 );// check in R: rnorm(2000,logit(0.1), 3) %>% inv_logit() %>% dens()
  
  // I_ini determined the season timing and certainly be a very low value
  logit(SIR_ini_mu[1,2]) ~ normal( logit(0.000002) , 0.7 ); // check in R: rnorm(2000,logit(0.000002),0.4) %>% inv_logit() %>% dens()
  logit(reciprocal_phi) ~ normal( logit(0.99) , 0.1 ); // check in R: rnorm(2000,logit(0.99),0.1) %>% inv_logit() %>% dens()
  
  // init priors
  logit(SIR_ini_mu[1,1]) ~ normal( logit(0.85) , 0.2 ); // check in R: rnorm(2000,logit(0.85),0.2) %>% inv_logit() %>% dens()
  logit(SIR_ini_mu[1,3]) ~ normal( logit(0.15) , 0.2 ); // check in R: rnorm(2000,logit(0.0015),0.4) %>% inv_logit() %>% dens()
  
  // priors on dispersion parameters
  sigma_prop_ili    ~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean LOW mean
  sigma_s              ~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean low mean
  sigma_i              ~ exponential(0.5); // parameter is the exponential RATE, so BIG numbers mean low mean
  
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
  array[n_week_fit] real<lower=0> delta_ili_abs_weekly_sum;
  array[n_week_fit, n_age_groups] int<lower=0> gen_ili_obs_fit;
  array[n_week_fit] int<lower=0> gen_ili_obs_fit_sum;
  // note: stan does not have 3-dimensional matrices, thus opting for arrays or 2-dimensional matrixes
  // note: matrix[n,m] M[o] creates an array of length o, each element contraining an nxm matrix, M CONFUSINGLY has then dimension [o,n,m]
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_ili_u[n_scenario];
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_ili_v[n_scenario];
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_ili_u_abs[n_scenario];
  matrix<lower=0>[n_day_project, n_age_groups] gen_delta_ili_v_abs[n_scenario];
  array[n_scenario, n_week_project, n_age_groups ] int<lower=0> gen_ili_u_obs_project; // unvaccinated
  array[n_scenario, n_week_project, n_age_groups ] int<lower=0> gen_ili_v_obs_project; // vaccinated
  array[n_scenario, n_week_project, n_age_groups ] int<lower=0> gen_ili_t_obs_project; // total
  array[n_scenario, n_week_project  ] int<lower=0> gen_ili_u_obs_project_sum;
  array[n_scenario, n_week_project  ] int<lower=0> gen_ili_v_obs_project_sum;
  array[n_scenario, n_week_project  ] int<lower=0> gen_ili_t_obs_project_sum;
  // gen_delta_ili_u_abs_daily
  array[n_scenario, n_day_project, n_age_groups] real gen_delta_ili_u_abs_daily; // unvaccinated
  array[n_scenario, n_day_project, n_age_groups] real gen_delta_ili_v_abs_daily; // vaccinated
  //
  array[n_scenario, n_week_project, n_age_groups] real gen_delta_ili_u_abs_weekly; // unvaccinated
  array[n_scenario, n_week_project, n_age_groups] real gen_delta_ili_v_abs_weekly; // vaccinated
  real Rnull_eff[n_season];
  real beta_noise;
  
  beta_noise = normal_rng( 1, 0.25 ); // to be applied as factor to log2( beta ), sd=1 interpreted as halfing / doubling of beta
  
  // --------------------------------simulate some quantities of interest
  // computation of Rnull_eff that is not quite correct due to age-structure
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  // --------------------------------simulate fitted observations
  // simulate fitted observations
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      gen_ili_obs_fit[t,a] = neg_binomial_2_rng( delta_ili_abs_weekly[t,a]+1e-6, phi );
    }
    gen_ili_obs_fit_sum[t] = sum( gen_ili_obs_fit[t,] );
    delta_ili_abs_weekly_sum[t] = sum( delta_ili_abs_weekly[t, ] );
  }
  
  // --------------------------------simulate projected observations
  for (j in 1:n_scenario) {
    // settings for the scenarios
    // define 2 local variables
    real beta_j; // scenario-specific beta
    matrix[n_day_project, n_age_groups] delta_vax_j; // scenario-specific vaccine uptake
    if ( axis_transmission[j]==0 ) beta_j = 1.0*beta; // more var: 1.0*2^( log2(beta)*beta_noise ) // pessimistic transmission
    if ( axis_transmission[j]==1 ) beta_j = 0.9*beta; // more var: 0.9*2^( log2(beta)*beta_noise ) // pessimistic transmission
    if ( axis_transmission[j]==2 ) beta_j = 1.1*beta; // more var: 1.1*2^( log2(beta)*beta_noise ) // optimisitc transmission
    if ( axis_vax[j]==0 ) delta_vax_j = delta_vax_real;
    if ( axis_vax[j]==1 ) delta_vax_j = delta_vax_opti;
    if ( axis_vax[j]==2 ) delta_vax_j = delta_vax_pess;
    if ( axis_vax[j]==3 ) delta_vax_j = delta_vax_null;
    int daily_counter;
    
    // Declare temporary variables outside the loop
    real gen_prev_S_u[n_age_groups];
    real gen_prev_I_u[n_age_groups];
    real gen_prev_R_u[n_age_groups];
    real gen_prev_S_v[n_age_groups];
    real gen_prev_I_v[n_age_groups];
    real gen_prev_R_v[n_age_groups];
    //
    real gen_curr_S_u[n_age_groups];
    real gen_curr_I_u[n_age_groups];
    real gen_curr_R_u[n_age_groups];
    real gen_curr_S_v[n_age_groups];
    real gen_curr_I_v[n_age_groups];
    real gen_curr_R_v[n_age_groups];
    //
    real curr_delta_ili_u[n_age_groups];
    real curr_delta_ili_v[n_age_groups];
    real curr_delta_ili_u_abs[n_age_groups];
    real curr_delta_ili_v_abs[n_age_groups];
    
    // main time loop
    daily_counter = 0;
    for (t in 1:n_multi_day_project) {
      // daily_counter: counting the number of days in relation to variable t
      if ((t-1)%n_daily_time_steps == 0){
        daily_counter = daily_counter+1;
      }
      //
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
        gen_curr_S_u[a] = SIR_ini[1,1,1];
        gen_curr_I_u[a] = SIR_ini[1,1,2];
        gen_curr_R_u[a] = 1 - (gen_curr_S_u[a] + gen_curr_I_u[a]);
        gen_curr_S_v[a] = 0;  // Adding a small number to avoid dividing by 0
        gen_curr_I_v[a] = 0;  // Adding a small number to avoid dividing by 0
        gen_curr_R_v[a] = 0;  // Adding a small number to avoid dividing by 0
        
        } else { // or update
        delta_infective_exposures_u = dt * beta_j * gen_prev_S_u[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(gen_prev_I_u)*1 + to_vector(gen_prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta_j * gen_prev_S_v[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(gen_prev_I_u)*1 + to_vector(gen_prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        //
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - gen_prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - gen_prev_I_v[a] * rate_infectious*dt;
        delta_R_u = gen_prev_I_u[a]*rate_infectious*dt;
        delta_R_v = gen_prev_I_v[a]*rate_infectious*dt;
        //
        gen_curr_S_u[a] = gen_prev_S_u[a] + delta_S_u - (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_prev_S_u[a] / (gen_prev_S_u[a] + gen_prev_R_u[a]);
        gen_curr_S_v[a] = gen_prev_S_v[a] + delta_S_v + (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_prev_S_u[a] / (gen_prev_S_u[a] + gen_prev_R_u[a]);
        gen_curr_I_u[a] = gen_prev_I_u[a] + delta_I_u;
        gen_curr_I_v[a] = gen_prev_I_v[a] + delta_I_v;
        gen_curr_R_u[a] = gen_prev_R_u[a] + delta_R_u - (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_prev_R_u[a] / (gen_prev_S_u[a] + gen_prev_R_u[a]);
        gen_curr_R_v[a] = gen_prev_R_v[a] + delta_R_v + (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_prev_R_u[a] / (gen_prev_S_u[a] + gen_prev_R_u[a]);
        //
        curr_delta_ili_u[a] = (delta_infective_exposures_u * 1) * prop_ili_mu[a];
        curr_delta_ili_v[a] = (delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili_mu[a];
        curr_delta_ili_u_abs[a] = curr_delta_ili_u[a] * pop_age_group[a,1];
        curr_delta_ili_v_abs[a] = curr_delta_ili_v[a] * pop_age_group[a,1];
        
        // Only store every n_daily_time_steps-th iteration
        if (t % n_daily_time_steps == 0) { 
          gen_S_u[t/n_daily_time_steps, a] = gen_curr_S_u[a];
          gen_S_v[t/n_daily_time_steps, a] = gen_curr_S_v[a];
          gen_I_u[t/n_daily_time_steps, a] = gen_curr_I_u[a];
          gen_I_v[t/n_daily_time_steps, a] = gen_curr_I_v[a];
          gen_R_u[t/n_daily_time_steps, a] = gen_curr_R_u[a];
          gen_R_v[t/n_daily_time_steps, a] = gen_curr_R_v[a];
          gen_delta_ili_u[j, t/n_daily_time_steps, a] = curr_delta_ili_u[a];
          gen_delta_ili_v[j, t/n_daily_time_steps, a] = curr_delta_ili_v[a];
          gen_delta_ili_u_abs[j, t/n_daily_time_steps, a] = curr_delta_ili_u_abs[a];
          gen_delta_ili_v_abs[j, t/n_daily_time_steps, a] = curr_delta_ili_v_abs[a];
        }
        
        if (t == 2*n_daily_time_steps) { // also impute the position of the first day
        gen_delta_ili_u_abs[j,1,a] = gen_delta_ili_u_abs[j,2,a];
        gen_delta_ili_v_abs[j,1,a] = gen_delta_ili_v_abs[j,2,a];
        gen_delta_ili_u[j,1,a] = gen_delta_ili_u[j,2,a];
        gen_delta_ili_v[j,1,a] = gen_delta_ili_v[j,2,a];
        }
        
        }
      } // n_age_groups loop
      // Update previous values with current values for next iteration
      gen_prev_S_u = gen_curr_S_u;
      gen_prev_I_u = gen_curr_I_u;
      gen_prev_R_u = gen_curr_R_u;
      gen_prev_S_v = gen_curr_S_v;
      gen_prev_I_v = gen_curr_I_v;
      gen_prev_R_v = gen_curr_R_v;
    } // n_day_project loop
  } // n_scenario loop
  
  // convert projections daily into weekly
  for (j in 1:n_scenario) {
    for (i in 1:n_week_project) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
        int day_end = day_start+6;
        gen_delta_ili_u_abs_weekly[j,i,a] = sum( gen_delta_ili_u_abs[j,day_start:day_end,a] );
        gen_delta_ili_v_abs_weekly[j,i,a] = sum( gen_delta_ili_v_abs[j,day_start:day_end,a] );
      }
    }
  }
  
  // simulate projections
  for (j in 1:n_scenario) {
    //for (j in 1:1) {
      for (t in 1:n_week_project) {
        for (a in 1:n_age_groups) {
          gen_ili_u_obs_project[j,t,a] = neg_binomial_2_rng( gen_delta_ili_u_abs_weekly[j,t,a]+1e-6 , phi ); // add small value to location parameter to avoid it being zero
          gen_ili_v_obs_project[j,t,a] = neg_binomial_2_rng( gen_delta_ili_v_abs_weekly[j,t,a]+1e-6 , phi ); // add small value to location parameter to avoid it being zero
          gen_ili_t_obs_project[j,t,a] = gen_ili_u_obs_project[j,t,a] + gen_ili_v_obs_project[j,t,a];
        }
        // sums across age-groups
        gen_ili_u_obs_project_sum[j,t] = sum(gen_ili_u_obs_project[j,t, ]);
        gen_ili_v_obs_project_sum[j,t]=  sum(gen_ili_v_obs_project[j,t, ]);
        gen_ili_t_obs_project_sum[j,t]=  sum(gen_ili_t_obs_project[j,t, ]);
      }
  }
  
  
}

