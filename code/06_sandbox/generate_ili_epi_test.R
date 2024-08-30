generate_ili_epi_test= function(par,stan_list){
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### DATA BLOCK  ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # data relevated for the fit 
  n_season = stan_list$n_season # number of seasons
  n_week_fit =  stan_list$n_week_fit # number of observable values, weekly
  n_day_fit = stan_list$n_day_fit # number of obervatble values, daily
  n_age_groups = stan_list$n_age_groups # number of age groups
  ili_obs_fit = stan_list$ili_obs_fit*0; # observed hospitalisations
  ili_obs_notna = stan_list$ili_obs_notna*0+1 # indicating non-missing data with 1, otherwise 0
  stan_list$ili_obs_notna = ili_obs_notna
  season_start = stan_list$season_start # indicating first week of a season with 1, the second week with 2, otherwise 0
  season_id_day = stan_list$season_id_day # indicating which seasn each obervable day belongs to
  pop = stan_list$pop # population size
  pop_age_group = stan_list$pop_age_group # population size per age group, requires to be a matrix 
  contact_matrix = stan_list$contact_matrix #contact matrix
  delta_vax = stan_list$delta_vax # daily fraction of newly vaccinated individuals per age group
  # epi parameters
  Rnull = stan_list$Rnull
  rate_infectious = stan_list$rate_infectious # infectious rate, such that beta = Rnull*rate_infectious
  # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf)
  ve_spread = stan_list$ve_spread # vaccine effectiveness on onward transmission/infectiousness
  ve_inf = stan_list$ve_inf # vaccine effectiveness on susceptability
  ve_ili_cond_inf = stan_list$ve_ili_cond_inf # vaccine effectiveness on severity, given infection
  # data relevant for projected scenarios
  n_week_project = stan_list$n_week_project # number of projected weeks
  n_day_project = stan_list$n_day_project # number of projected days
  n_scenario = stan_list$n_scenario # number of projected scenarios
  axis_transmission = stan_list$axis_transmission # indicator for the transmission scenario axis
  axis_vax = stan_list$axis_vax # indicator for the vaccine scenario axis
  delta_vax_real = stan_list$delta_vax_real # daily assumed vax uptake in projection period
  delta_vax_opti = stan_list$delta_vax_opti # daily assumed vax uptake in projection period
  delta_vax_pess = stan_list$delta_vax_pess # daily assumed vax uptake in projection period
  delta_vax_null = stan_list$delta_vax_null # daily assumed vax uptake in projection period
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### TRANSFORMED DATA BLOCK  ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  beta = rate_infectious * Rnull
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### PARAMETER BLOCK ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SIR_ini_mu = c(0.77, 0.000003 , 1 - 0.77 - 0.000003) # overall season mean 
  SIR_ini = t(matrix(c(
    c(0.77, 0.000003, NA), 
    c(0.77, 0.000003, NA),
    c(0.77, 0.000003, NA)), nrow=3)) # S I R initial values per season, 1 can be replaced by n_age_groups
  SIR_ini[,3] = 1 - (SIR_ini[,1] + SIR_ini[,2] )
  #
  prop_ili_mu = c(0.10)
  prop_ili = t(matrix(rep(0.10,n_age_groups*n_season), nrow=n_age_groups))
  # Hyper parameter structure
  set.seed(seed=params$simulation_seed)
  sigma_prop_ili_age = 2
  sigma_prop_ili_season = 0.5
  sigma_i_season = 2
  if (F) {
    prop_ili_age_factor = rnorm(n=n_age_groups,0,sd=sigma_prop_ili_age) %>% exp()
    prop_ili_season_factor = rnorm(n=n_season,0,sd=sigma_prop_ili_season) %>% exp()
    for (age_i in 1:n_age_groups) {
      for (season_i in 1:n_season) {
        prop_ili[season_i, age_i] = prop_ili_mu * prop_ili_season_factor[season_i] * prop_ili_age_factor[age_i] 
      }
    }
    #
    i_season_factor = rnorm(n=n_season,0,sd=sigma_i_season) %>% exp()
    for (season_i in 1:n_season) SIR_ini[season_i,2] = SIR_ini_mu[2] * i_season_factor[season_i]
    SIR_ini[,3] = 1 - (SIR_ini[,1] + SIR_ini[,2] ) # the R compartment absorbs the change in I
  }
  #
  reciprocal_phi = 0.05 # overdipersion parameter for ili obs fit, var=mu+reciprocal_phi*mu^2
  sim_par = 
    list(
      SIR_ini_mu=SIR_ini_mu,
      SIR_ini = SIR_ini,
      prop_ili_mu = prop_ili_mu,
      prop_ili = prop_ili,
      sigma_prop_ili_age = sigma_prop_ili_age,
      sigma_prop_ili_season = sigma_prop_ili_season,
      sigma_i_season = sigma_i_season,
      reciprocal_phi = reciprocal_phi
      
    )
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### TRANSFORMED PARAMETER BLOCK  ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
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
        S_u[t,a] = SIR_ini[season_id_day[t], 1] # * pop_age_group[a, 1] / pop # rescaling 
        I_u[t,a] = SIR_ini[season_id_day[t], 2] # * pop_age_group[a, 1] / pop # rescaling
        R_u[t,a] = SIR_ini[season_id_day[t], 3] # * pop_age_group[a, 1] / pop # rescaling
        S_v[t,a] = 0  # at start of season, no one is vaccinated
        I_v[t,a] = 0  # at start of season, no one is vaccinated
        R_v[t,a] = 0  # at start of season, no one is vaccinated
      }
      
    } else {
      for(a in 1:n_age_groups){  
        delta_infective_exposures_u = beta * S_u[t-1,a]  * sum( contact_matrix[a ,] * ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_spread)) )
        delta_infective_exposures_v = beta * S_v[t-1,a]  * sum( contact_matrix[a ,] * ( I_u[t-1,]*1 + I_v[t-1,]*(1-ve_spread)) ) * (1 - ve_inf)
        
        delta_S_u = -delta_infective_exposures_u
        delta_S_v = -delta_infective_exposures_v
        delta_I_u = delta_infective_exposures_u - I_u[t-1,a] * rate_infectious
        delta_I_v = delta_infective_exposures_v - I_v[t-1,a] * rate_infectious 
        delta_R_u = I_u[t-1,a]*rate_infectious
        delta_R_v = I_v[t-1,a]*rate_infectious 
        
        # infection
        S_u[t,a] = S_u[t-1,a] + delta_S_u 
        S_v[t,a] = S_v[t-1,a] + delta_S_v 
        I_u[t,a] = I_u[t-1,a] + delta_I_u
        I_v[t,a] = I_v[t-1,a] + delta_I_v
        R_u[t,a] = R_u[t-1,a] + delta_R_u
        R_v[t,a] = R_v[t-1,a] + delta_R_v 
        # vaccination
        S_u[t,a] = S_u[t,a] - data.frame(delta_vax)[t,a] * S_u[t,a] 
        S_v[t,a] = S_v[t,a] + data.frame(delta_vax)[t,a] * S_u[t,a] 
        R_u[t,a] = R_u[t,a] - data.frame(delta_vax)[t,a] * R_u[t,a] 
        R_v[t,a] = R_v[t,a] + data.frame(delta_vax)[t,a] * R_u[t,a] 
        
        delta_ili[t,a] = (delta_infective_exposures_u * 1 + delta_infective_exposures_v * (1-ve_ili_cond_inf) ) * prop_ili[season_id_day[t], a]
        delta_ili_abs[t,a] = delta_ili[t,a] * pop_age_group[a,1]
        
        if (season_start[t]==2){
          # fill first position of the season in other vectors
          delta_ili_abs[t-1,a] = delta_ili_abs[t,a]
          delta_ili[t-1,a] = delta_ili[t,a]
        }
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
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### LIKELIHOOD/PRIOR BLOCK  ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  set.seed(seed= params$simulation_seed)
  # nonoise = T
  for (t in 1:n_week_fit) {
    for (a in 1:n_age_groups) {
      if (ili_obs_notna[t,a]==1) ili_obs_fit[t,a] = rnbinom(1, mu=delta_ili_abs_weekly[t,a]+1e-6, size=phi )
      # if (nonoise) if (ili_obs_notna[t,a]==1) ili_obs_fit[t,a] = round(delta_ili_abs_weekly[t,a])
      }
  }
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### GENERATED QUANTITIES BLOCK   ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  stan_list$sim_par = sim_par
  for (age_group_i in 1:n_age_groups){
    stan_list$ili_obs_fit[,age_group_i] = ili_obs_fit[,age_group_i]
  }
  # summary stats
  stan_list$cum_ili_obs_log = rowsum(x=stan_list$ili_obs_fit,group=stan_list$season_id_week,na.rm = T) %>% rowSums() %>% log()
  stan_list$n_ili_obs_notna = rowsum(x=stan_list$ili_obs_notna,group=stan_list$season_id_week,na.rm = T) %>% rowSums()
  
  return(stan_list)
}

