if (F){
  df_Rt_without_imm = compute_Rt_without_immunity(par=res_optim$par,
                                                  Rt_seasonal,
                                                  for_optim = F,
                                                  p=p,
                                                  days_to_week=days_to_week_sim,
                                                  sim_o = list(
                                                    inf_ini = 1000*4,
                                                    Imm_1_ini = Imm_1_ini,
                                                    Imm_2_ini = Imm_2_ini,
                                                    si_distribution = p$si_discr_dist
                                                  ) )
  
  
  df_Rt_without_imm %>% g()
  df_Rt_without_imm$Rt_seasonal %>% plot()
  df_Rt_without_imm$Rt_eff %>% plot()
  df_Rt_without_imm$Imm %>% plot()
  
  df_Rt_without_imm %>% ggplot(aes(x=date)) +
    geom_line(aes(y=value))
}


compute_Rt_without_immunity = function( par,
                                        # ability to fix parameters
                                        wane_imm_fract_weekly_1_log=NA,
                                        Imm_2_partial=NA,
                                        wane_imm_fract_weekly_2_log=NA,
                                        #
                                        df_Rt_value_pop,
                                        for_optim=F,
                                        p,
                                        days_to_week,
                                        sim_o=NA) {
  # df_Rt_value_pop must be daily
  # we apply a chunking-of-days technique to speed up optim etc
  simulate = !is.na(sim_o)[1] # logical to determine whether or not to forward simulate 
  ## recasting parameters (needed due to optim use)
  p$ar_and_re = exp(par[1])
  wane_imm_fract_weekly_1_logit = wane_imm_fract_weekly_1_log
  if (is.na(wane_imm_fract_weekly_1_log)) wane_imm_fract_weekly_1_logit = par[2]
  Imm_2_partial_logit = Imm_2_partial
  if (is.na(Imm_2_partial)) Imm_2_partial_logit = par[3]
  wane_imm_fract_weekly_2_logit = wane_imm_fract_weekly_2_log
  if (is.na(wane_imm_fract_weekly_2_log)) wane_imm_fract_weekly_2_logit = par[4]
  
  
  # waning compartments could be further extended
  # wane_imm_fract_weekly_3_logit = par[4]
  # Imm_3_partial
  df_all  = df_Rt_value_pop
  df_all = df_all[1:days_to_week$n_full_weeks,]
  # weeks and days
  nweeks = days_to_week$nweeks
  weekly_daygroups_t = days_to_week$weekly_daygroups_t
  weekly_daygroups_tlast = days_to_week$weekly_daygroups_tlast
  
  ## transformations
  Imm_2_partial = inv_logit(Imm_2_partial_logit)
  wane_imm_fract_weekly_1 = inv_logit(wane_imm_fract_weekly_1_logit) # scale daily/weekly
  wane_imm_fract_weekly_2 = inv_logit(wane_imm_fract_weekly_2_logit) # scale daily/weekly
  # wane_imm_fract_weekly_3 = exp(wane_imm_fract_weekly_3_log)
  
  ## country loop
  # for each location, simulate immunity and update target based on conditions
  location_v = unique(df_all$location_name)
  df_save = list()
  target = 0 # initiate target for optimisation
  min_target = 0 # initiate minimum target for optimisation
  for (i in location_v) { # i = location_v[1]
    df_i = df_all %>% filter(location_name==i) 
    
    #### value: infection intensity
    if (!simulate) {
      # normalise infection intensity (so that infection attack rate is 1% each year)
      laste_year_ii = mlast(df_i$value,n=365)/mlast(df_i$population,n=365)
      df_i$value = df_i$value * ( 0.01/sum(  laste_year_ii  )  )
    }
    if (simulate) {
      df_i$value = 0; df_i$value[ weekly_daygroups_tlast[[1]][1] ] = sim_o$inf_ini
    }
    
    #### Imm_X: the different immunity compartments
    df_i$Imm_1 = NA; 
    df_i$Imm_2 = NA; 
    if (!simulate) {
      # initial Rt and immunity
      Rt_ini = mean(df_i$Rt[1:p$Rt_for_ini_est]) # avoiding impact of noise in Rt 
      Imm_ini = (p$R0 - Rt_ini)/p$R0 # back-calculating immunity from Rt and R0
      # assuming that all immunity is recent (thus all into first immunity compartment)
      df_i$Imm_1[weekly_daygroups_tlast[[1]]] = Imm_ini # full immunity
      df_i$Imm_2[weekly_daygroups_tlast[[1]]] = 0 # partial immunity
    }
    if (simulate) {
      # initial immunity
      df_i$Imm_1[weekly_daygroups_tlast[[1]]] = sim_o$Imm_1_ini
      df_i$Imm_2[weekly_daygroups_tlast[[1]]] = sim_o$Imm_2_ini
    }
    # build immunity from the immunity compartments
    df_i$Imm = NA; 
    df_i$Imm[weekly_daygroups_tlast[[1]]] = df_i$Imm_1[weekly_daygroups_tlast[[1]]] + df_i$Imm_2[weekly_daygroups_tlast[[1]]]*Imm_2_partial # immunity combined
    
    #### waning_X: waning from the different immunity compartments
    # initiate waning
    df_i$waning_1 = NA; df_i$waning_1[weekly_daygroups_tlast[[1]]] = 0
    df_i$waning_2 = NA; df_i$waning_2[weekly_daygroups_tlast[[1]]] = 0
    df_i$new_Imm = NA; df_i$new_Imm[weekly_daygroups_tlast[[1]]] = 0
    
    #### Rt_without_immunity and Rt with immmunity
    if (!simulate){
      # compute Rt without immunity
      df_i$Rt_without_immunity = NA; 
      df_i$Rt_without_immunity[weekly_daygroups_tlast[[1]]] = df_i$Rt[weekly_daygroups_tlast[[1]]] * 1/(1-Imm_ini)
    }
    if (simulate) {
      # Rt without immunity exist, compute Rt effective
      df_i$Rt_eff = NA
      df_i$Rt_without_immunity = df_i$Rt_seasonal 
      df_i$Rt_eff[weekly_daygroups_tlast[[1]]] = 
        df_i$Rt_without_immunity[weekly_daygroups_tlast[[1]]] * (1-df_i$Imm[weekly_daygroups_tlast[[1]]])
    }
    
    #### sX: different scenarios
    ## counter-factuals/scenarios/helpers
    df_i$s1 = NA; df_i$s1[weekly_daygroups_tlast[[1]]] = 1
    df_i$s1_1 = NA; df_i$s1_1[weekly_daygroups_tlast[[1]]] = 1
    df_i$s1_2 = NA; df_i$s1_2[weekly_daygroups_tlast[[1]]] = 0
    # 
    df_i$s2 = NA; df_i$s2[weekly_daygroups_tlast[[1]]] = 1
    #
    df_i$s3 = NA; df_i$s3[weekly_daygroups_tlast[[1]]] = df_i$s1[weekly_daygroups_tlast[[1]]]*df_i$Imm[weekly_daygroups_tlast[[1]]]
    
    # time loop (through chunks of days of a given location)
    for (t in 1:nweeks ) { # t = 1
      # deal with weekly to daily
      t_days  = weekly_daygroups_t[[t]]
      tlast_days  = weekly_daygroups_tlast[[t]]
      # simlate infections (value) daily stemming from last week 
      if (simulate) {
        # handle daily growth simulation 1 day-chunk at a time
        for (chunk_day_i in 1:length(tlast_days) ){ # chunk_day_i = 1
          whichday = tlast_days[chunk_day_i]
          current_infection = df_i$value[whichday]
          current_Rt = df_i$Rt_eff[whichday]
          current_offspring = current_infection * current_Rt
          
          si_v_along = c(1:length(sim_o$si_distribution))
          
          ## add offspring incident infections back into value
          distributed_offspring = current_offspring * sim_o$si_distribution
          # truncate if offspring goes beyond 
          n_too_long = ( whichday+mlast(si_v_along) ) - days_to_week$n_full_weeks
          trunk_v = c( 1:( length(whichday+si_v_along)-max(n_too_long,0) ) ) 
          # add the offspring
          df_i$value[(whichday+si_v_along)[trunk_v]] = df_i$value[(whichday+si_v_along)[trunk_v]] + 
            distributed_offspring[trunk_v]
        } # df_i$value[1: (whichday+mlast(si_v_along)) ] %>% plot()
      }
      
      # immunity waning (allowing different immunity compartments)
      # Imm_1 is the full immunity compartment
      # Imm_2 3, ... are partial immunity compartments
      waning_1 = df_i$Imm_1[tlast_days[1]]*wane_imm_fract_weekly_1 
      waning_2 = df_i$Imm_2[tlast_days[1]]*wane_imm_fract_weekly_2 
      df_i$waning_1[t_days] = waning_1/7 # distribute over days
      df_i$waning_2[t_days] = waning_2/7 # distribute over days
      # accumulation of new immunity (proportional to: value = "infection intensity",dark_factor,NE)
      new_infections = sum(df_i$value[tlast_days]/df_i$population[tlast_days])
      new_Imm = new_infections * p$ar_and_re
      df_i$new_Imm[t_days] = new_Imm/7 # distribute over days
      # new infections also come from partly susceptible Imm compartments
      full_susc = 1 - ( df_i$Imm_1[tlast_days[1]] + df_i$Imm_2[tlast_days[1]] )
      if (full_susc<0&for_optim==F) {print("warning0")}
      new_Imm_from_Imm_2 = df_i$Imm_2[tlast_days[1]]*Imm_2_partial / 
        ( df_i$Imm_2[tlast_days[1]]*Imm_2_partial + full_susc*1 )
      # apply changes to Immunity compartments
      df_i$Imm_1[t_days] = df_i$Imm_1[tlast_days[1]] - waning_1 + new_Imm
      df_i$Imm_2[t_days] = df_i$Imm_2[tlast_days[1]] - waning_2 + waning_1 - new_Imm_from_Imm_2*new_Imm
      if (any(df_i$Imm_1[t_days]<0)&for_optim==F) print("warning1")
      if (any(df_i$Imm_2[t_days]<0)&for_optim==F) {print("warning2")}
      # recompute combined "effective" immunity
      df_i$Imm[t_days] = df_i$Imm_1[t_days] + df_i$Imm_2[t_days]*Imm_2_partial
      # compute new Rt without immunity (Rt AS IF immunity would have vanished)
      if (!simulate) df_i$Rt_without_immunity[t_days] = df_i$Rt[t_days] * 1/(1-df_i$Imm[t_days])
      if ( simulate) df_i$Rt_eff[t_days] = df_i$Rt_without_immunity[t_days] * (1-df_i$Imm[t_days])
      ## scenarios 
      # 1: just waning shape
      waning_s1_1 = df_i$s1_1[tlast_days[1]]*wane_imm_fract_weekly_1 
      waning_s1_2 = df_i$s1_2[tlast_days[1]]*wane_imm_fract_weekly_2
      df_i$s1_1[t_days] = df_i$s1_1[tlast_days] - waning_s1_1 
      df_i$s1_2[t_days] = df_i$s1_2[tlast_days] - waning_s1_2 + waning_s1_1
      df_i$s1[t_days] = df_i$s1_1[t_days] + df_i$s1_2[t_days]*Imm_2_partial
      # 2: just accounting for new immunity
      new_week_eff_immunity = sum(df_i$value[tlast_days]/df_i$population[tlast_days]) * p$ar_and_re
      new_week_Rt_modifier = 1 - new_week_eff_immunity
      new_week_Rt_modifier_inv = 1 / new_week_Rt_modifier
      df_i$s2[t_days] = df_i$s2[tlast_days] * new_week_Rt_modifier_inv
      # 3: just accounting for immunity from the start
      df_i$s3[t_days] = df_i$s1[t_days]*df_i$Imm[1]
    }
    
    if (!simulate){
      ## "target" allows multiple constraints
      # regularize: immunity last year ~ immunity 1 year later
      Imm_diff = df_i$Imm[df_i$Rt_at_one==1] - df_i$Imm[df_i$Rt_at_one==2]
      target_Imm_diff = -dnorm(x=Imm_diff,mean = 0,sd=0.05, log = T)
      target_Imm_diff_min = -dnorm(x=0,mean = 0,sd=0.05, log = T)
      # regularize: R0 based on Rt and R0
      # note: Rt_without immunity can exceed R0, since Rt and Immunity beyong the initial condition don't know about each other 
      max_Rt_without_immunity = max(df_i$Rt_without_immunity,na.rm=T)
      target_Rt_without_immunity = -dnorm(x=max_Rt_without_immunity,mean = p$R0,sd=1, log = T)
      target_Rt_without_immunity_min = -dnorm(x=p$R0,mean = p$R0,sd=1, log = T)
      # regularize; effective waning 
      
      # could also regularize; max immunity in the population
      
      target = target + target_Imm_diff + target_Rt_without_immunity
      min_target = min_target = target_Imm_diff_min + target_Rt_without_immunity_min # -1.157855
    }
    
    # print(paste("min target:",min_target %>% round(2)))
    # print(paste("Max new_Imm",max(df_i$new_Imm) %>% round(3),i))
    # print(paste("Average waning",mean(df_i$waning) %>% round(3),i))
    df_save[[i]] = df_i
  }
  # pack it all up
  df_return = bind_rows(df_save) %>% 
    mutate(target=target, 
           min_target=min_target)
  
  # return depending on where the function is used
  if (for_optim==T) return(target)
  if (for_optim==F) return(df_return)
}
