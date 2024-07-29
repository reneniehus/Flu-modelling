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
  // test simplex[3] SIR_ini[n_season, 1]; // S I R initial values per season, 1 can be replaced by n_age_groups
  simplex[3] SIR_ini_scaled[n_season, 1]; // S I R initial values per season, 1 can be replaced by n_age_groups
  simplex[3] SIR_ini_mu[1]; // overall season mean 
  // test real<lower=0, upper=1> prop_ili[n_season, n_age_groups]; // proportion of infections that are ili 
  real<lower=0, upper=1> prop_ili_scaled[n_season, n_age_groups]; // proportion of infections that are ili 
  real<lower=0, upper=1> prop_ili_mu[n_age_groups]; // overall mean over season 
  // dispersion parameters
  real<lower=0> sigma_prop_ili;
  real<lower=0> sigma_s;
  real<lower=0> sigma_i;
  
  real<lower=0, upper=1> reciprocal_phi; // overdipersion parameter for ili obs fit, var=mu+reciprocal_phi*mu^2
}

transformed parameters {
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
  // array[n_day_fit,n_age_groups] real<lower=0> delta_ili_abs_daily; // ili/detectable incidence in absolute numbers per day
  real phi; // dispersion parameter, var=mu+reciprocal_phi*mu^2
  // weekly stuff
  array[n_week_fit,n_age_groups] real<lower=0> delta_ili_abs_weekly; // ili/detectable incidence in absolute numbers, weekly aggregate
  
  // Overdispersion
  phi = 1 / reciprocal_phi; // dispersion parameter: var=mu+reciprocal_phi*mu^2
  
  // test
  real<lower=0, upper=1> prop_ili[n_season, n_age_groups];
  array[n_season,1,3] real<lower=0, upper=1> SIR_ini;
  real sum_val;
  for (i in 1:n_season) {
    for (j in 1:n_age_groups) {
      prop_ili[i, j] = prop_ili_scaled[i, j] * 0.2;
    }
    SIR_ini[i,1,1] = SIR_ini_scaled[i,1,1];
    SIR_ini[i,1,2] = SIR_ini_scaled[i,1,2] * 0.01;
    SIR_ini[i,1,3] = SIR_ini_scaled[i,1,3];
    sum_val = SIR_ini[i,1,1] + SIR_ini[i,1,2] + SIR_ini[i,1,3];
    SIR_ini[i,1,1] = SIR_ini[i,1,1] / sum_val;
    SIR_ini[i,1,2] = SIR_ini[i,1,2] / sum_val;
    SIR_ini[i,1,3] = SIR_ini[i,1,3] / sum_val;
  }
  
  
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
    // print(t);
    // print(n_multi_day_fit);
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
    // print("10");
    
    // If it is first day of season AND first iteration of the day
    if ( season_start[daily_counter_fit[t]]==1 && ((t-1)%n_daily_time_steps)==0 ){
      // print("0a");
      //
      // initiate the compartments based on current season\
      // S I R initial values age dist corrected
      for(a in 1:n_age_groups){
        // print("0b");
        curr_S_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 1] * pop_age_group[a, 1] / pop; // rescaling
        curr_I_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 2] * pop_age_group[a, 1] / pop; // rescaling
        curr_R_u[a] = SIR_ini[season_id[daily_counter_fit[t]], 1, 3] * pop_age_group[a, 1] / pop; // rescaling
        curr_S_v[a] = 0;  // at start of season, no one is vaccinated
        curr_I_v[a] = 0;
        curr_R_v[a] = 0;
        curr_delta_ili[a] = 1;
        curr_delta_ili_abs[a] = 1;
        // print("Current values");
        // print(t);
        // print(curr_S_u);
        // print(curr_I_u);
        // print(curr_R_u);
        // print("====");
        
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
      // print("0c");
      
    } else {
      for(a in 1:n_age_groups){
        // print("----");
        // print(a);
        delta_infective_exposures_u = dt * beta * prev_S_u[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta * prev_S_v[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        
        // print("1+++++++");
        // print(t);
        // print(n_daily_time_steps);
        // print(1.0/n_daily_time_steps);
        // print("----");
        // print(dt);
        // print(beta);
        // print(prev_S_u[a]);
        // print("----");
        // print(delta_infective_exposures_u);
        // print(delta_infective_exposures_v);
        // print("----");
        // print(sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ));
        // print(to_vector(contact_matrix[ : , a]'));
        // print((to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        // print("----");
        // print(to_vector(contact_matrix[ : , a]'));
        // print(to_vector(prev_I_u));
        // print(to_vector(prev_I_v));
        // print((to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        // print("1-------");
        
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        delta_I_u = delta_infective_exposures_u - prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - prev_I_v[a] * rate_infectious*dt;
        delta_R_u = prev_I_u[a]*rate_infectious*dt;
        delta_R_v = prev_I_v[a]*rate_infectious*dt;
        // print("2+++++++++");
        // print(delta_S_u);
        // print(delta_S_v);
        // print("HERE:");
        // print(delta_I_u);
        // print(delta_I_v);
        // print(delta_R_u);
        // print(delta_R_v);
        // print("-");
        // print(rate_infectious);
        // print(prev_I_u[a]);
        // print("-");
        // print((delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps));
        // print("2-------");
        
        
        curr_S_u[a] = prev_S_u[a] + delta_S_u - (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_S_v[a] = prev_S_v[a] + delta_S_v + (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_I_u[a] = prev_I_u[a] + delta_I_u;
        curr_I_v[a] = prev_I_v[a] + delta_I_v;
        curr_R_u[a] = prev_R_u[a] + delta_R_u - (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_R_v[a] = prev_R_v[a] + delta_R_v + (delta_vax[daily_counter_fit[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        
        // print("3+++++++");
        // print(delta_infective_exposures_u);
        // print(delta_infective_exposures_v);
        // print(prop_ili[season_id[daily_counter_fit[t]], a] );
        // print("3-------");
        curr_delta_ili[a] = (delta_infective_exposures_u * 1 + delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili[season_id[daily_counter_fit[t]], a];
        curr_delta_ili_abs[a] = curr_delta_ili[a] * pop_age_group[a,1];
        
        // if (t==2){
          //   reject("t==2, done");
          // }
          // if (max(curr_I_u) > 0.01) {
            //   print("curr_I_u:");
            //   print(t);
            //   print(curr_I_u);
            // }
            // if (max(curr_I_u) > 0.10) {
              //   print("curr_I_u is large enough - 3!");
              // }
              //
              // if (max(curr_I_u) > 0.05) {
                //   print("curr_I_u is large enough - 2!");
                // }
                
                if (t % n_daily_time_steps == 0) { // Only store every n_daily_time_steps-th iteration
                // print("+++++++++");
                // print(t);
                // print(a);
                // print(curr_S_u);
                // print(curr_S_v);
                // print(curr_I_u);
                // print(curr_I_v);
                // print(curr_R_u);
                // print(curr_R_v);
                // print("--");
                // print(curr_delta_ili);
                // print(curr_delta_ili_abs);
                // print("-------------");
                // if (max(curr_I_u) > 0.01) {
                  //   reject("curr_I_u is large enough?");
                  // }
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
  
  
  
  /*
  print("5a");
  for (i in 1:n_multi_day_fit) {
  for (j in 1:n_age_groups) {
  if (is_nan(delta_ili[i,j]) || delta_ili[i,j] < 0) {
  print(i);
  print(j);
  reject("delta_ili[", i, ",", j, "] is invalid: ", delta_ili[i,j]);
  //reject("rejecting");
  }
  }
  }
  */
  
  //
  // print("5");
  // print(delta_ili_abs);
  // print("5a");
  /*
  // convert projections multi-daily into daily
  for (i in 1:n_day_fit) {
  for (a in 1:n_age_groups) {
  // define 2 local variables
  int day_start = (i-1)*n_daily_time_steps+1; // f(i=1)=1 , f(i=2)=8
  int day_end = day_start+(n_daily_time_steps-1);
  delta_ili_abs_daily[i,a] = sum( delta_ili_abs[day_start:day_end,a] );
  delta_ili_abs_daily[i,a] = sum( delta_ili_abs[day_start:day_end,a] );
  }
  }
  */
  //print("6");
  //print(delta_ili_abs);
  //print("6a");
  
  // convert daily to weekly
  for (i in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      // define 2 local variables
      int day_start = (i-1)*7+1;
      int day_end = day_start+6;
      delta_ili_abs_weekly[i,a] = sum( delta_ili_abs[day_start:day_end,a] );
    }
  }
  //print("7");
  //print(delta_ili_abs_weekly);
  
}
//print("8");
//}

model {
  
  // --------------------------------likelihood part
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      // print("running neg_binomial_2");
      // print(t);
      // print(a);
      // print(delta_ili_abs_weekly[t,a]);
      // print("end neg_binomial_2");
      if (ili_obs_notna[t,a]==1) ili_obs_fit[t,a] ~ neg_binomial_2( delta_ili_abs_weekly[t,a]+1e-6, phi ) ; // TODO: remove this line and see if you get priors back
    }
  }
  
  // --------------------------------prior part
  // prior
  // real alpha14;
  // real alpha24;
  // real mu4;
  // real s24;
  // mu4 = 0.015;
  // s24 = 0.5*mu4*(1-mu4);
  // alpha14 = mu4 * (mu4 * (1 - mu4) / s24 - 1);
  // alpha24 = (1 - mu4) * (mu4 * (1 - mu4) / s24 - 1);
  for (a in 1:n_age_groups) {
    //logit( prop_ili[,a] ) ~ normal( logit( prop_ili_mu[a] ) , sigma_prop_ili );
    //prop_ili[,a]   ~ beta( alpha14, alpha24 );
    prop_ili_scaled[,a] ~ uniform(0, 1);
  }
  //logit( SIR_ini[,1,1] )   ~ normal( logit( SIR_ini_mu[1,1] )   , sigma_s );
  // logit( SIR_ini[,1,2] )   ~ normal( logit( SIR_ini_mu[1,2] )   , sigma_i );
  //logit( SIR_ini[,1,3] )   ~ normal( logit( SIR_ini_mu[1,3] )   , sigma_s );
  
  // for (a in 1:n_age_groups) {
    //   prop_ili[,a] ~ normal( logit( prop_ili_mu[a] ) , sigma_prop_ili );
    // }
    
    // Calculate alpha(=alpha1) and beta(=alpha2) from mean (mu) and variance (s2)
    // real alpha11;
    // real alpha21;
    // real mu1;
    // real s21;
    // mu1 = 0.85;
    // s21 = 0.1*mu1*(1-mu1);
    // alpha11 = mu1 * (mu1 * (1 - mu1) / s21 - 1);
    // alpha21 = (1 - mu1) * (mu1 * (1 - mu1) / s21 - 1);
    // SIR_ini[,1,1]   ~ beta( alpha11, alpha21 );
    
    // real alpha12;
    // real alpha22;
    // real mu2;
    // real s22;
    // mu2 = 0.000002;
    // s22 = 0.1*mu2*(1-mu2);
    // alpha12 = mu2 * (mu2 * (1 - mu2) / s22 - 1);
    // alpha22 = (1 - mu2) * (mu2 * (1 - mu2) / s22 - 1);
    //SIR_ini[,1,2]   ~ beta( alpha12, alpha22 );
    
    // Uniform priors for SIR_ini
    //SIR_ini[,1,1]   ~ uniform( 0.8, 0.9 );
    SIR_ini_scaled[,1,2]   ~ uniform( 0, 1 );
    //SIR_ini[,1,2]   ~ uniform( 0.0000005, 0.000005 );
    //vector[3] alpha = to_vector([1, 1e4, 1e3]);
    //SIR_ini[,1,] ~ dirichlet(alpha);
    
    // real alpha13;
    // real alpha23;
    // real mu3;
    // real s23;
    // mu3 = 0.15;
    // s23 = 0.1*mu3*(1-mu3);
    // alpha13 = mu3 * (mu3 * (1 - mu3) / s23 - 1);
    // alpha23 = (1 - mu3) * (mu3 * (1 - mu3) / s23 - 1);
    // SIR_ini[,1,3]   ~ beta( alpha13, alpha23 );
    
    // SIR_ini[,1,1]   ~ beta(  10  , 1.75 );
    // //SIR_ini[,1,2]   ~ beta( 0.000001999992   , 0.999994 );
    // SIR_ini[,1,2]   ~ beta( alpha1, alpha2 );
    // SIR_ini[,1,3]   ~ beta( 2   , 10.5 );
    
    //
    // Z2 ~ normal(0,1);
    // logit( SIR_ini[,1,2] ) ~ logit(SIR_ini_mu[1,2]) + Z2 * sigma_i; // non-centered parametrization
    // // SIR_ini_logit[:,2] = logit(SIR_ini_mu[1,2]) + Z2 * sigma_s; // Apply non-centered parametrization
    // // SIR_ini[:,2] ~ normal( inv_logit(SIR_ini_logit[:,2]), 1); // Transform back to the original scale using the inverse logit function
    // Z3 ~ normal(0,1);
    // logit( SIR_ini[,1,3] )   = Z3 * sigma_s + logit(SIR_ini_mu[1,3]); // non-centered parametrization
    //
    
    // more priors: put priors on things that we want to fixate more
    // logit(prop_ili_mu) ~ normal( logit(0.1) , 3 );// check in R: rnorm(2000,logit(0.1), 3) %>% inv_logit() %>% dens()
    logit(prop_ili_mu) ~ normal( logit(0.015) , 0.5 );// check in R: rnorm(2000,logit(0.1), 3) %>% inv_logit() %>% dens()
    
    // I_ini determined the season timing and certainly be a very low value
    // R and S we want to keep free to allow learning about immunity
    // logit(SIR_ini_mu[1,2]) ~ normal( logit(0.0015) , 0.4 ); // check in R: rnorm(2000,logit(0.0015),0.4) %>% inv_logit() %>% dens()
    // logit(reciprocal_phi) ~ normal( logit(0.99) , 0.1 ); // check in R: rnorm(2000,logit(0.99),0.1) %>% inv_logit() %>% dens()
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
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_S_u;
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_I_u;
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_R_u;
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_S_v;
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_I_v;
  matrix<lower=0,upper=1>[n_multi_day_project, n_age_groups] gen_R_v;
  array[n_week_fit] real<lower=0> delta_ili_abs_weekly_sum;
  array[n_week_fit, n_age_groups] int<lower=0> gen_ili_obs_fit;
  array[n_week_fit] int<lower=0> gen_ili_obs_fit_sum;
  // note: stan does not have 3-dimensional matrices, thus opting for arrays or 2-dimensional matrixes
  // note: matrix[n,m] M[o] creates an array of length o, each element contraining an nxm matrix, M CONFUSINGLY has then dimension [o,n,m]
  matrix<lower=0>[n_multi_day_project, n_age_groups] gen_delta_ili_u[n_scenario];
  matrix<lower=0>[n_multi_day_project, n_age_groups] gen_delta_ili_v[n_scenario];
  matrix<lower=0>[n_multi_day_project, n_age_groups] gen_delta_ili_u_abs[n_scenario];
  matrix<lower=0>[n_multi_day_project, n_age_groups] gen_delta_ili_v_abs[n_scenario];
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
  
  // print("dims(gen_delta_ili): ", dims(gen_delta_ili) ); //
  // print("dims(gen_delta_ili_abs): ", dims(gen_delta_ili_abs) ); // [6,364,1]
  // print("dims(gen_ili_obs_project): ", dims(gen_ili_obs_project) ); //
  // print("dims(gen_delta_ili_abs_weekly): ", dims(gen_delta_ili_abs_weekly) ); //
  
  // --------------------------------simulate some quantities of interest
  // computation of Rnull_eff that is not quite correct due to age-structure
  for (season_i in 1:n_season) {
    Rnull_eff[season_i] = Rnull*(1-sum(SIR_ini[ season_i,,3 ]));
  }
  
  // --------------------------------simulate fitted observations
  // simulate fitted observations
  //print("simulate fitted observations");
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {

      gen_ili_obs_fit[t,a] = neg_binomial_2_rng( delta_ili_abs_weekly[t,a]+1e-6, phi );
      //print("simulate fitted observations 2");
    }
    gen_ili_obs_fit_sum[t] = sum( gen_ili_obs_fit[t,] );
    delta_ili_abs_weekly_sum[t] = sum( delta_ili_abs_weekly[t, ] );
    //print("simulate fitted observations 3");
  }
  // print("end simulate fitted observations");
  
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
    real prev_S_u[n_age_groups];
    real prev_I_u[n_age_groups];
    real prev_R_u[n_age_groups];
    real prev_S_v[n_age_groups];
    real prev_I_v[n_age_groups];
    real prev_R_v[n_age_groups];
    //
    real curr_S_u[n_age_groups];
    real curr_I_u[n_age_groups];
    real curr_R_u[n_age_groups];
    real curr_S_v[n_age_groups];
    real curr_I_v[n_age_groups];
    real curr_R_v[n_age_groups];
    //
    real curr_delta_ili_u[n_age_groups];
    real curr_delta_ili_v[n_age_groups];
    real curr_delta_ili_u_abs[n_age_groups];
    real curr_delta_ili_v_abs[n_age_groups];
    
    // main time loop
    //print("main time loop");
    //print(j);
    //print("above is scenario j");
    daily_counter = 0;
    for (t in 1:n_multi_day_project) {
      
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
        
        //print("main time loop t==1");
        curr_S_u[a] = SIR_ini_mu[1,1];
        curr_I_u[a] = SIR_ini_mu[1,2];
        curr_R_u[a] = 1 - (curr_S_u[t,a] + curr_I_u[t,a]);
        curr_S_v[a] = 0;  // Adding a small number to avoid dividing by 0
        curr_I_v[a] = 0;  // Adding a small number to avoid dividing by 0
        curr_R_v[a] = 0;  // Adding a small number to avoid dividing by 0
        
        // old time variables - can delete
        // //print("main time loop t==1");
        // gen_S_u[t,a] = SIR_ini_mu[1,1];
        // gen_I_u[t,a] = SIR_ini_mu[1,2];
        // gen_R_u[t,a] = 1 - (gen_S_u[t,a] + gen_I_u[t,a]);
        // gen_S_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        // gen_I_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        // gen_R_v[t,a] = 0;  // Adding a small number to avoid dividing by 0
        
        } else { // or update
        //
        //print("main time loop 2");
        // delta_infective_exposures_u = dt * beta_j * gen_S_u[t-1,a] *sum(contact_matrix[ : , a]' .* (gen_I_u[t-1,] + gen_I_v[t-1,]*(1-ve_spread)));
        // delta_infective_exposures_v = dt * beta_j * gen_S_v[t-1,a] *sum(contact_matrix[ : , a]' .* (gen_I_u[t-1,] + gen_I_v[t-1,]*(1-ve_spread))) * (1-ve_inf);
        delta_infective_exposures_u = dt * beta_j * prev_S_u[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) );
        delta_infective_exposures_v = dt * beta_j * prev_S_v[a] * sum(to_vector(contact_matrix[ : , a]') .* (to_vector(prev_I_u)*1 + to_vector(prev_I_v)*(1-ve_spread)) ) * (1 - ve_inf);
        
        
        delta_S_u = -delta_infective_exposures_u;
        delta_S_v = -delta_infective_exposures_v;
        //print("aa");
        // delta_I_u = delta_infective_exposures_u - gen_I_u[t-1,a] * rate_infectious*dt;
        // delta_I_v = delta_infective_exposures_v - gen_I_v[t-1,a] * rate_infectious*dt;
        delta_I_u = delta_infective_exposures_u - prev_I_u[a] * rate_infectious*dt;
        delta_I_v = delta_infective_exposures_v - prev_I_v[a] * rate_infectious*dt;
        
        //print("aa2");
        // delta_R_u = gen_I_u[t-1,a] * rate_infectious*dt;
        // delta_R_v = gen_I_v[t-1,a] * rate_infectious*dt;
        delta_R_u = prev_I_u[a]*rate_infectious*dt;
        delta_R_v = prev_I_v[a]*rate_infectious*dt;
        //
        //print("aa3");
        //print(t);
        //print(daily_counter-1);
        //print(a);
        //print(delta_vax_j[daily_counter-1,a]);
        //print(gen_S_u[t-1,a]);
        //print(gen_R_u[t-1,a]);
        //print("aa3a");
        curr_S_u[a] = prev_S_u[a] + delta_S_u - (delta_vax_j[daily_counter[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_S_v[a] = prev_S_v[a] + delta_S_v + (delta_vax_j[daily_counter[t-1],a]/n_daily_time_steps) * prev_S_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_I_u[a] = prev_I_u[a] + delta_I_u;
        curr_I_v[a] = prev_I_v[a] + delta_I_v;
        curr_R_u[a] = prev_R_u[a] + delta_R_u - (delta_vax_j[daily_counter[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        curr_R_v[a] = prev_R_v[a] + delta_R_v + (delta_vax_j[daily_counter[t-1],a]/n_daily_time_steps) * prev_R_u[a] / (prev_S_u[a] + prev_R_u[a]);
        
        
        // gen_S_u[t,a] = gen_S_u[t-1,a] + delta_S_u - (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        // gen_S_v[t,a] = gen_S_v[t-1,a] + delta_S_v + (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_S_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        //print("aa4");
        // gen_I_u[t,a] = gen_I_u[t-1,a] + delta_I_u;
        // gen_I_v[t,a] = gen_I_v[t-1,a] + delta_I_v;
        //print("aa5");
        // gen_R_u[t,a] = gen_R_u[t-1,a] + delta_R_u - (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        // gen_R_v[t,a] = gen_R_v[t-1,a] + delta_R_v + (delta_vax_j[daily_counter,a]/n_daily_time_steps) * gen_R_u[t-1,a] / (gen_S_u[t-1,a] + gen_R_u[t-1,a]);
        //
        //print("aa6");
        curr_delta_ili_u[a] = (delta_infective_exposures_u * 1) * prop_ili_mu[a];
        curr_delta_ili_v[a] = (delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili_mu[a];
        curr_delta_ili_u_abs[a] = curr_delta_ili_u[a] * pop_age_group[a,1];
        curr_delta_ili_v_abs[a] = curr_delta_ili_v[a] * pop_age_group[a,1];
        
        // gen_delta_ili_u[j,t,a] = (delta_infective_exposures_u * (1-0) ) * prop_ili_mu[a];
        // gen_delta_ili_v[j,t,a] = (delta_infective_exposures_v * (1-ve_ili_cond_inf)) * prop_ili_mu[a];
        // 
        // gen_delta_ili_u_abs[ j,t,a] = gen_delta_ili_u[j,t,a] * pop_age_group[a,1] ;
        // gen_delta_ili_v_abs[ j,t,a] = gen_delta_ili_v[j,t,a] * pop_age_group[a,1] ;
        //
        //print("aa7");
        
        if (t % n_daily_time_steps == 0) { // Only store every n_daily_time_steps-th iteration
        gen_S_u[t/n_daily_time_steps, a] = curr_S_u[a];
        gen_S_v[t/n_daily_time_steps, a] = curr_S_v[a];
        gen_I_u[t/n_daily_time_steps, a] = curr_I_u[a];
        gen_I_v[t/n_daily_time_steps, a] = curr_I_v[a];
        gen_R_u[t/n_daily_time_steps, a] = curr_R_u[a];
        gen_R_v[t/n_daily_time_steps, a] = curr_R_v[a];
        gen_delta_ili_u[n_scenario, t/n_daily_time_steps, a] = curr_delta_ili_u[a];
        gen_delta_ili_v[n_scenario, t/n_daily_time_steps, a] = curr_delta_ili_v[a];
        gen_delta_ili_u_abs[n_scenario, t/n_daily_time_steps, a] = curr_delta_ili_u_abs[a];
        gen_delta_ili_v_abs[n_scenario, t/n_daily_time_steps, a] = curr_delta_ili_v_abs[a];
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
      prev_S_u = curr_S_u;
      prev_I_u = curr_I_u;
      prev_R_u = curr_R_u;
      prev_S_v = curr_S_v;
      prev_I_v = curr_I_v;
      prev_R_v = curr_R_v;
    } // n_day_project loop
  } // n_scenario loop
  
  // convert projections multi-daily into daily
  for (j in 1:n_scenario) {
    for (i in 1:n_day_project) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (i-1)*n_daily_time_steps+1; // f(i=1)=1 , f(i=2)=8
        int day_end = day_start+(n_daily_time_steps-1);
        gen_delta_ili_u_abs_daily[j,i,a] = sum( gen_delta_ili_u_abs[j,day_start:day_end,a] );
        gen_delta_ili_v_abs_daily[j,i,a] = sum( gen_delta_ili_v_abs[j,day_start:day_end,a] );
      }
    }
  }
  
  // convert projections daily into weekly
  for (j in 1:n_scenario) {
    for (i in 1:n_week_project) {
      for (a in 1:n_age_groups) {
        // define 2 local variables
        int day_start = (i-1)*7+1; // f(i=1)=1 , f(i=2)=8
        int day_end = day_start+6;
        gen_delta_ili_u_abs_weekly[j,i,a] = sum( gen_delta_ili_u_abs_daily[j,day_start:day_end,a] );
        gen_delta_ili_v_abs_weekly[j,i,a] = sum( gen_delta_ili_v_abs_daily[j,day_start:day_end,a] );
      }
    }
  }
  
  // simulate projections
  for (j in 1:n_scenario) {
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

