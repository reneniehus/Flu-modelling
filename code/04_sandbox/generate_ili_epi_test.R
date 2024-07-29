generate_ili_epi_test= function(stan_list){
  # daily stuff
  # SIR compartments unvaccinated and vaccinated
  # variables starting with delta_ are incidence variables that are daily, unless it is specificed otherwise (e.g. through _weekly )
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_u; // susceptible compartment, relative to population size, unvaccinated
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_u; // infetious compartment,   relative to population size, unvaccinated
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_u; // recovered compartment,   relative to population size, unvaccinated
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] S_v; // susceptible compartment, relative to population size, vaccinated
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] I_v; // infetious compartment,   relative to population size, vaccinated
  #matrix<lower=0, upper=1>[n_day_fit,n_age_groups] R_v; // recovered compartment,   relative to population size, vaccinated
  #array[n_day_fit,n_age_groups] real<lower=0, upper=1> delta_ili; // ili/detectable incidence relative to population size
  #array[n_day_fit,n_age_groups] real<lower=0> delta_ili_abs; // ili/detectable incidence in absolute numbers
  #real phi; // dispersion parameter, var=mu+reciprocal_phi*mu^2
  #// weekly stuff
  #array[n_week_fit,n_age_groups] real<lower=0> delta_ili_abs_weekly; // ili/detectable incidence in absolute numbers, weekly aggregate
  
  # Data
  # data relevated for the fit 
  n_season = stan_list$n_season # number of seasons
  n_week_fit =  stan_list$n_week_fit # number of observable values, weekly
  n_day_fit = stan_list$n_day_fit # number of obervatble values, daily
  n_age_groups = stan_list$n_age_groups # number of age groups
  ####ili_obs_fit[n_week_fit, n_age_groups]; // observed hospitalisations
  ili_obs_notna = stan_list$ili_obs_notna # indicating non-missing data with 1, otherwise 0
  season_start = stan_list$season_start # indicating first week of a season with 1, the second week with 2, otherwise 0
  season_id = stan_list$season_id # indicating which seasn each obervable day belongs to
  pop = stan_list$pop # population size
  pop_age_group = stan_list$pop_age_group # population size per age group, requires to be a matrix 
  contact_matrix = stan_list$contact_matrix #contact matrix
  delta_vax = stan_list$delta_vax # daily fraction of newly vaccinated individuals per age group
  # data relevant for projected scenarios
  n_week_project = stan_list$n_week_project # number of projected weeks
  n_day_project = stan_list$n_day_project # number of projected days
  #n_scenario = stan_list$n_scenario # number of projected scenarios
  #axis_transmission = stan_list$axis_transmission # indicator for the transmission scenario axis
  #axis_vax = stan_list$axis_vax # indicator for the vaccine scenario axis
  #delta_vax_real = stan_list$delta_vax_real # daily assumed vax uptake in projection period
  #delta_vax_opti = stan_list$delta_vax_opti # daily assumed vax uptake in projection period
  #delta_vax_pess = stan_list$delta_vax_pess # daily assumed vax uptake in projection period
  #delta_vax_null = stan_list$delta_vax_null # daily assumed vax uptake in projection period
  
  # epi parameters
  Rnull = stan_list$Rnull
  rate_infectious = stan_list$rate_infectious # infectious rate, such that beta = Rnull*rate_infectious
  # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf)
  ve_spread = stan_list$ve_spread # vaccine effectiveness on onward transmission/infectiousness
  ve_inf = stan_list$ve_inf # vaccine effectiveness on susceptability
  ve_ili_cond_inf = stan_list$ve_ili_cond_inf # vaccine effectiveness on severity, given infection
  # epi parameters
  beta = rate_infectious * Rnull
  
  
  # parameters to be 'fitted'
  # note: a simplex of length 3, has 2 free parameters
  # note: simplex[n] X[m,o] creates an m x o sized array of simplex, each of size n
  SIR_ini = t(matrix(c(c(0.85, 0.000001, 0.149999), 
              c(0.84, 0.000003, 0.159997),
              c(0.87, 0.000002, 0.129998)), nrow=3)) # S I R initial values per season, 1 can be replaced by n_age_groups
  #SIR_ini_mu = c(0.78, 0.02, 0.2) # overall season mean 
  prop_ili = t(matrix(rep(0.01,16), nrow=4))
  prop_ili_mu = rep(0.01,4)
  #prop_ili = t(matrix(rep(c(0.05, 0.01, 0.002, 0.002),4), nrow=4)) # proportion of infections that are ili 
  #prop_ili_mu = c(0.05, 0.01, 0.002, 0.002) # overall mean over season 

  reciprocal_phi = 0.98 # overdipersion parameter for ili obs fit, var=mu+reciprocal_phi*mu^2
  
  # Overdispersion
  phi = 1 / reciprocal_phi # dispersion parameter: var=mu+reciprocal_phi*mu^2
  
  # Initialize few varibles
  S_u = matrix(NA, n_day_fit, n_age_groups)
  I_u = matrix(NA, n_day_fit, n_age_groups)
  R_u = matrix(NA, n_day_fit, n_age_groups)
  S_v = matrix(NA, n_day_fit, n_age_groups)
  I_v = matrix(NA, n_day_fit, n_age_groups)
  R_v = matrix(NA, n_day_fit, n_age_groups)
  delta_ili = matrix(NA, n_day_fit, n_age_groups)
  delta_ili_abs = matrix(NA, n_day_fit, n_age_groups)
  delta_ili_abs_weekly = matrix(NA, n_week_fit, n_age_groups)
  
  # loop through all days
  for (t in 1:n_day_fit){
    
    if ( season_start[t]==1 ){
      # initiate the compartments based on current season\
      # S I R initial values age dist corrected
      for(a in 1:n_age_groups){
        S_u[t,a] = SIR_ini[season_id[t], 1] * pop_age_group[a, 1] / pop # rescaling 
        I_u[t,a] = SIR_ini[season_id[t], 2] * pop_age_group[a, 1] / pop # rescaling
        R_u[t,a] = SIR_ini[season_id[t], 3] * pop_age_group[a, 1] / pop # rescaling
        S_v[t,a] = 0  # at start of season, no one is vaccinated
        I_v[t,a] = 0 
        R_v[t,a] = 0 
        
      }
      
    } else {
      for(a in 1:n_age_groups){  
        delta_infective_exposures_u = beta * S_u[t-1,a]  * sum(contact_matrix[ , a] * ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_spread)) )
        delta_infective_exposures_v = beta * S_v[t-1,a]  * sum(contact_matrix[ , a] * ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_spread)) ) * (1 - ve_inf)
        
        delta_S_u = -delta_infective_exposures_u
        delta_S_v = -delta_infective_exposures_v
        delta_I_u = delta_infective_exposures_u - I_u[t-1,a] * rate_infectious
        delta_I_v = delta_infective_exposures_v - I_v[t-1,a] * rate_infectious 
        delta_R_u = I_u[t-1,a]*rate_infectious
        delta_R_v = I_v[t-1,a]*rate_infectious 
        
        S_u[t,a] = S_u[t-1,a] + delta_S_u - data.frame(delta_vax)[t-1,a] * S_u[t-1,a] #/ (S_u[t-1,a] + R_u[t-1,a])
        S_v[t,a] = S_v[t-1,a] + delta_S_v + data.frame(delta_vax)[t-1,a] * S_u[t-1,a] #/ (S_u[t-1,a] + R_u[t-1,a])
        I_u[t,a] = I_u[t-1,a] + delta_I_u
        I_v[t,a] = I_v[t-1,a] + delta_I_v
        R_u[t,a] = R_u[t-1,a] + delta_R_u - data.frame(delta_vax)[t-1,a] * R_u[t-1,a] #/ (S_u[t-1,a] + R_u[t-1,a])
        R_v[t,a] = R_v[t-1,a] + delta_R_v + data.frame(delta_vax)[t-1,a] * R_u[t-1,a] #/ (S_u[t-1,a] + R_u[t-1,a])
        
        #
        delta_ili[t,a] = (delta_infective_exposures_u * 1 + delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili[season_id[t], a]
        delta_ili_abs[t,a] = delta_ili[t,a] * pop_age_group[a,1]
        
        if (season_start[t]==2){
          # fill first position of the season in other vectors
          delta_ili_abs[t-1,a] = delta_ili_abs[t,a]
          delta_ili[t-1,a] = delta_ili[t,a]
        }
      }
      
      
    } # end of daily loop
    
    # convert daily to weekly
    for (i in 1:n_week_fit) {
      for (a in 1:n_age_groups) {
        # define 2 local variables
        day_start = (i-1)*7+1; 
        day_end = day_start+6;
        delta_ili_abs_weekly[i,a] = sum( delta_ili_abs[day_start:day_end,a] );
      }
    }
    
  }
  
  delta_ili_abs_weekly = round(delta_ili_abs_weekly, digits = 2)
  # Output
  stan_list$ili_obs_fit$age_00_04 = delta_ili_abs_weekly[,1] %>% as.integer()
  stan_list$ili_obs_fit$age_05_14 = delta_ili_abs_weekly[,2] %>% as.integer()
  stan_list$ili_obs_fit$age_15_64 = delta_ili_abs_weekly[,3] %>% as.integer()
  stan_list$ili_obs_fit$age_65_99 = delta_ili_abs_weekly[,4] %>% as.integer()
  #stan_list$ili_obs_fit$age_total = rowSums(delta_ili_abs_weekly) %>% as.integer()

  return(stan_list)
  
}

