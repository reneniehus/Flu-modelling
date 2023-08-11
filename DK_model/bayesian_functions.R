wave_settings_SIR = list(
  I_ini=0.0001,
  S_ini= 0.8,
  beta = 0.20,
  prop_severe = 0.02,
  sigma_m = log(50)*0.0001
)
sim_wave_SIR = function(wave_settings_SIR,sim_name="no_name"){
  
  # simulate 1 wave of a pathogen
  # strongly simplifying assumptions: well mixed
  pop=1e6 # population size, big enough to make things look smooth
  t_v = c(1:500) # 
  na_v = t_v*NA
  severe_obs = na_v
  # epi parameters
  beta = wave_settings_SIR$beta; logit(beta)
  prop_severe = wave_settings_SIR$prop_severe; logit(prop_severe)
  I_ini = wave_settings_SIR$I_ini 
  S_ini = wave_settings_SIR$S_ini 
  sigma_m = wave_settings_SIR$sigma_m # 0.01% noise
  rate_infectious = 0.1042 # 1/(7+2.6)
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.1 wave
  # start with fully susceptible population
  # introduce hypothetical Omicron BA.1 (at fraction I_ini into population) and then save the resulting level of immunity
  
  I[1] = I_ini
  S[1] = S_ini # 
  R[1] = 1 - (I_ini+S_ini)
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("Wave not complete!")
  }
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures,sim_name) -> wave_df
  if (F){
    mcaption = paste("Contagion wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    wave_df %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") + labs(caption=mcaption)
  }
  
  
  return(wave_df)
}

wave_settings_SIR_Rnull = list(
  I_ini=0.0001,
  S_ini= 0.8,
  Rnull = 2.5,
  prop_severe = 0.02,
  sigma_m = log(50)*0.01
)
sim_wave_SIR_Rnull = function(wave_settings_SIR_Rnull,sim_name="no_name"){
  
  # simulate 1 wave of a pathogen
  # strongly simplifying assumptions: well mixed
  pop=1e6 # population size, big enough to make things look smooth
  t_v = c(1:500) # 
  na_v = t_v*NA
  severe_obs = na_v
  # epi parameters
  Rnull = wave_settings_SIR_Rnull$Rnull
  rate_infectious = 0.1042 # 1/(7+2.6)
  beta = rate_infectious*Rnull; logit(beta)
  prop_severe = wave_settings_SIR_Rnull$prop_severe; logit(prop_severe)
  I_ini = wave_settings_SIR_Rnull$I_ini ; logit(I_ini)
  S_ini = wave_settings_SIR_Rnull$S_ini ; logit(S_ini)
  sigma_m = wave_settings_SIR_Rnull$sigma_m # 0.01% noise
  
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.1 wave
  # start with fully susceptible population
  # introduce hypothetical Omicron BA.1 (at fraction I_ini into population) and then save the resulting level of immunity
  
  I[1] = I_ini
  S[1] = S_ini # 
  R[1] = 1 - (I_ini+S_ini)
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta * S[t-1] * I[t-1]                     # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1] * rate_infectious      # eq3
    delta_R[t] = I[t-1] * rate_infectious                                     # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t] * pop                                    # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("Wave not complete!")
  }
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures,sim_name) -> wave_df
  if (F){
    mcaption = paste("Contagion wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    wave_df %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") + labs(caption=mcaption)
  }
  
  
  return(wave_df)
}

sim_wave_SIR_Rnull_w = function(wave_settings_SIR_Rnull,sim_name="no_name"){
  
  # simulate 1 wave of a pathogen
  # strongly simplifying assumptions: well mixed
  pop=1e6 # population size, big enough to make things look smooth
  n_week = 30
  n_day = n_week * 7
  t_v = c(1:n_day) # 
  t_vw = c(1:n_week)
  na_v = t_v*NA
  na_vw = t_vw*NA
  severe_obs = na_v
  # epi parameters
  Rnull = wave_settings_SIR_Rnull$Rnull
  rate_infectious = 0.2777778 # 1/(3.6)
  beta = rate_infectious*Rnull; logit(beta)
  prop_severe = wave_settings_SIR_Rnull$prop_severe; logit(prop_severe)
  I_ini = wave_settings_SIR_Rnull$I_ini ; logit(I_ini)
  S_ini = wave_settings_SIR_Rnull$S_ini ; logit(S_ini)
  sigma_m = wave_settings_SIR_Rnull$sigma_m # 0.01% noise
  
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  severe_mean_weekly = na_vw
  severe_mean_weekly_fat = na_vw
  severe_obs_weekly = na_vw
  severe_obs_weekly_fat = na_vw
  i_day = na_vw
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.1 wave
  # start with fully susceptible population
  # introduce hypothetical Omicron BA.1 (at fraction I_ini into population) and then save the resulting level of immunity
  
  I[1] = I_ini
  S[1] = S_ini # 
  R[1] = 1 - (I_ini+S_ini)
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("Wave not complete!")
  }
  # fill first position in other vectors
  severe_obs[1] = severe_obs[2]
  severe_mean[1] = severe_mean[2]
  severe_obs_log[1] = log(severe_obs[1])
  
  # make wave fatter
  severe_mean %>% sum() -> mtarget
  severe_fat = severe_mean + 2 * c(rep(0,1*14),severe_mean[ 1: (length(severe_mean)-1*14) ]) + 
    1 * c(rep(0,2*14),severe_mean[ 1: (length(severe_mean)-2*14) ])
  severe_fat = severe_fat/sum(severe_fat)
  severe_fat = severe_fat*mtarget
  
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_fat,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures,sim_name) -> wave_df_d
  
  # convert daily to weekly
  for (i in 1:n_week) {
    day_start = (i-1)*7+1 # f(i=1)=1 , f(i=2)=8
    day_end = day_start+6
    severe_mean_weekly[i] = sum( severe_mean[day_start:day_end] )
    severe_mean_weekly_fat[i] = sum( severe_fat[day_start:day_end] )
    severe_obs_weekly[i] = rpois(n=1,lambda=severe_mean_weekly[i])
    severe_obs_weekly_fat[i] = rpois(n=1,lambda=severe_mean_weekly_fat[i])
    i_day[i] = t_v[day_start]
  }
  tibble(t_vw,i_day,severe_mean_weekly,severe_obs_weekly,severe_mean_weekly_fat,severe_obs_weekly_fat,pop) -> wave_df_w
  
  if (F){
    mcaption = paste("Contagion wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    wave_df_d %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") +
      geom_point(data=wave_df_w,aes(x=i_day,y=severe_obs_weekly/7),col="red")+
      labs(caption=mcaption)
  }
  
  
  return(wave_df_w)
}

wave_settings = list(
  pop=1e6,
  n_week = 30,
  I_ini=0.0001,
  S_ini= 0.8,
  Rnull = 2.5,
  rate_infectious = 0.2777778, # 1/(3.6),
  prop_severe = 0.02,
  seed_rand = 12,
  add_observation_settings=list(sigma_m=log(50)*0.01,
                                seed=12)
)
# helper functions for sim_wave_SIR_ll
make_observation_from_epi_signal = function(severe_mean, add_observation_settings){
  set.seed(add_observation_settings$seed)
  severe_obs_log = rnorm( length(severe_mean),mean = log( severe_mean ),sd = add_observation_settings$sigma_m )
  severe_obs = exp(severe_obs_log)
}

sim_wave_SIR_ll = function(wave_settings,sim_name="no_name"){
  # simulate 1 wave of a pathogen
  # strongly simplifying assumptions: well mixed
  pop=wave_settings$pop # population size, big enough to make things look smooth
  n_week = wave_settings$n_week
  n_day = n_week * 7
  t_v = c(1:n_day) # 
  t_vw = c(1:n_week)
  na_v = t_v*NA
  na_vw = t_vw*NA
  severe_obs = na_v
  # epi parameters
  Rnull = wave_settings$Rnull
  rate_infectious = wave_settings$rate_infectious
  beta = rate_infectious*Rnull; logit(beta)
  prop_severe = wave_settings$prop_severe; logit(prop_severe)
  I_ini = wave_settings$I_ini ; logit(I_ini)
  S_ini = wave_settings$S_ini ; logit(S_ini)
  sigma_m = wave_settings$sigma_m # 0.01% noise
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  severe_mean_weekly = na_vw
  severe_mean_weekly_fat = na_vw
  severe_obs_weekly = na_vw
  severe_obs_weekly_fat = na_vw
  i_day = na_vw
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  I[1] = I_ini
  S[1] = S_ini 
  R[1] = 1 - (I_ini+S_ini)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("Wave not complete!")
  }
  
  # make wave fatter (overlapping sub-epidemics)
  if (wave_settings$fat_wave==T) {
    # have 1 wave first, then 2 waves, then 1 wave again, then normalise
    severe_mean %>% sum() -> mtarget
    severe_fat = severe_mean + 2 * c(rep(0,1*14),severe_mean[ 1: (length(severe_mean)-1*14) ]) + 
      1 * c(rep(0,2*14),severe_mean[ 1: (length(severe_mean)-2*14) ])
    severe_fat = severe_fat/sum(severe_fat)
    severe_fat = severe_fat*mtarget
    severe_mean = severe_fat
  }
  
  # make observation process based on epi signal
  id_second_to_last = (2:length(severe_mean))
  severe_obs[id_second_to_last] = make_observation_from_epi_signal( severe_mean[id_second_to_last], wave_settings$add_observation_settings )
  # fill first position in other vectors
  # FIXME: maybe these should not be imputed, but simply not fit!
  severe_obs[1] = severe_obs[2]
  severe_mean[1] = severe_mean[2]
  
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures,sim_name) -> wave_df_d
  
  
  # convert time series to line list data
  df_for_ll = wave_df_d %>% select(t_v,severe_obs)
  wave_ll = make_ll_from_epi( df_for_ll , make_ll_settings )
  
  make_ll_from_epi = function(df_for_ll , make_ll_settings) {
    #
    #
  }
  
  df_out = list(
    wave_df_daily=wave_df_d,
    wave_ll = wave_ll
    
  )
  
  
  if (F){
    mcaption = paste("Contagion wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    wave_df_d %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") +
      geom_point(data=wave_df_w,aes(x=i_day,y=severe_obs_weekly/7),col="red")+
      labs(caption=mcaption)
  }
  
  return(df_out)
}


sim_wave_DK = function(wave_settings){
  
  # BA.1 wave parameters
  # since the BA.1 wave is not fitted to epi data, the prop_severe and noise of observation system (sigma_m) is NOT fitted, it does not impact the final immunity
  # only beta impacts the final immunity
  # during the time where BA.1 and BA.2 waves overlap - we only fit BA.1 
  pop=1e6 # population size
  t_v = c(1:200) # 
  na_v = t_v*NA
  severe_obs = na_v
  # epi parameters
  beta = 0.15; logit(beta)
  prop_severe = 0.02; logit(prop_severe)
  I_ini = 0.001 
  sigma_m = log(50)*0.0001 # 0.01% noise
  rate_infectious = 0.1042 # 1/(7+2.6)
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.1 wave
  # start with fully susceptible population
  # introduce hypothetical Omicron BA.1 (at fraction I_ini into population) and then save the resulting level of immunity
  
  I[1] = I_ini
  S[1] = 1 - I_ini # 
  R[1] = 0 
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("BA.1 wave not complete!")
  }
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures) -> ba1_hosp
  if (F){
    mcaption = paste("BA.1 wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    ba1_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") + labs(caption=mcaption)
  }
  
  # BA.2 wave parameters - we consider all immunity with respect to BA.2!
  # thus the "artificial" BA.1 wave is created to produce sufficient immunity against BA.2
  # making the BA.1 specific beta also an "artificial" beta
  # happens in the same population
  # this time the "dark factor" (prop_severe) and the variability (sigma_m) of the observation process matters
  # the I_ini needs fitting
  # assume the same transmission parameter for this variant as for the BA.5
  pop=1e6 # population size
  t_v = c(1:200) # 
  na_v = t_v*NA
  severe_obs = na_v
  # epi parameters
  beta = 0.7; logit(beta)
  prop_severe = 0.02; logit(prop_severe)
  I_ini = 0.04 
  prop_immune_escape = 0 # proportion of recovered that becomes available due to immunune evasion
  sigma_m = log(50)*0.0001 # 0.01% noise
  rate_infectious = 0.1042 # 1/(7+2.6)
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.2 wave
  # population has susceptibility of the final stage of BA.1 wave
  S[1] = mlast(ba1_hosp$S) + mlast(ba1_hosp$R) * prop_immune_escape
  I[1] = I_ini
  R[1] = 1 - (S[1] + I[1]) 
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("BA.2 wave not complete!")
  }
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures) -> ba2_hosp
  if (F){
    mcaption = paste("BA.2 wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    ba2_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") + labs(caption=mcaption)
  }
  
  # BA.5 wave parameters
  # happens just like the BA.2 wave
  pop=1e6 # population size
  t_v = c(1:200) # 
  na_v = t_v*NA
  severe_obs = na_v
  # epi parameters
  # beta = 0.3; logit(beta)
  prop_severe = 0.02; logit(prop_severe)
  I_ini = 0.04 
  prop_immune_escape = 0.7 # proportion of recovered that becomes available due to immunune evasion
  sigma_m = log(50)*0.0001 # 0.01% noise
  rate_infectious = 0.1042 # 1/(7+2.6)
  
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions for BA.5 wave
  # population has susceptibility of the final stage of BA.1 wave
  S[1] = mlast(ba2_hosp$S) + mlast(ba2_hosp$R) * prop_immune_escape
  I[1] = I_ini
  R[1] = 1 - (S[1] + I[1]) 
  set.seed(12)
  for ( t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]                        # eq1
    delta_S[t] = -delta_infective_exposures[t]                                # eq2
    delta_I[t] = delta_infective_exposures[t] - I[t-1]*rate_infectious        # eq3
    delta_R[t] = I[t-1]*rate_infectious                                       # eq4
    S[t] = S[t-1] + delta_S[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq5
    severe_mean[t] = delta_severe[t]*pop                                      # eq6
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq7
    severe_obs[t] = exp(severe_obs_log[t])
    # warning if first wave is not complete 
    if (t==length(t_v) & I[t]>I_ini ) print("BA.5 wave not complete!")
  }
  # first wave saving and plotting
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,R,delta_S,delta_I,delta_R,delta_infective_exposures) -> ba5_hosp
  if (F){
    mcaption = paste("BA.5 wave, starting with", I_ini, "infections and simulated over" ,length(t_v), "days")
    ba2_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=severe_mean),col="blue") + labs(caption=mcaption)
  }
  bind_rows(
    ba1_hosp %>% mutate(voc="ba1"),
    ba2_hosp %>% mutate(voc="ba2"),
    ba5_hosp %>% mutate(voc="ba5")
  ) -> df_hosp
  
  return(df_hosp)
}

if (F){
  df_hosp = sim_wave_DK(NULL)
  df_hosp %>% filter(voc=="ba1") %>% 
    ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
    geom_point(alpha=0.2)+
    geom_line(aes(y=severe_mean),col="orange")+labs(subtitle = "BA.1 wave") -> p1 # orange because this epi curve is not meant to be realistic
  df_hosp %>% filter(voc=="ba2") %>% 
    ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
    geom_point(alpha=0.2)+
    geom_line(aes(y=severe_mean),col="blue") +labs(subtitle = "BA.2 wave") -> p2
  df_hosp %>% filter(voc=="ba5") %>% 
    ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
    geom_point(alpha=0.2)+
    geom_line(aes(y=severe_mean),col="blue") +labs(subtitle = "BA.5 wave") -> p3
  (p1 + p2 + p3) & coord_cartesian(ylim=c(0,800) )
}



sim_wave_SEIR = function(wave_settings){
  pop=1e6 # population size
  t_v = c(1:121) # 
  na_v = t_v*NA
  severe_obs = na_v
  severe_obs[2] = 50 # after dividing by pop: 0.00002
  
  beta = 0.3; logit(beta)
  prop_severe = 0.02; logit(prop_severe)
  S_ini = 0.80; logit(S_ini)
  sigma_m = log(50)*0.0001 # 0.01% noise
  rate_latent = 0.385
  rate_infectious = 0.142
  # initialising empty vectors
  delta_severe = na_v
  S = na_v
  I = na_v
  E = na_v
  R = na_v
  delta_infective_exposures = na_v
  severe_mean = na_v
  severe_obs_log = na_v
  delta_S = na_v
  delta_E = na_v
  delta_I = na_v
  delta_R = na_v
  
  # Initial conditions
  S[1] = S_ini # immunity at t1
  # Setting the initial I is a back-calculation from observed severe under some assumptions
  # * no noise in measurement
  severe_obs_log[2] = log(severe_obs[2])
  severe_mean[2] = severe_obs[2] # this is not exactly true
  delta_severe[2] = severe_mean[2]/pop
  delta_infective_exposures[2] = delta_severe[2]/prop_severe # due to eq6
  I[1] = delta_infective_exposures[2] / (beta*S[1] )         # due to eq1
  # E depends on the recent history of the SEIR model
  
  # let's assume steady-state for now
  E[1] = I[1]*rate_infectious/rate_latent # 2.6 %
  
  R[1] = 1 - (S[1] + I[1] + E[1]); # 30 %
  delta_severe[1] = severe_obs[1]/pop
  set.seed(12)
  for (t in 2:length(t_v) ) { # t = 2
    delta_infective_exposures[t] = beta*S[t-1] *I[t-1]             # eq1
    delta_S[t] = -delta_infective_exposures[t]                     # eq2
    delta_E[t] = delta_infective_exposures[t] - E[t-1]*rate_latent # eq3
    delta_I[t] = E[t-1]*rate_latent - I[t-1]*rate_infectious       # eq4
    delta_R[t] = I[t-1]*rate_infectious                            # eq5
    S[t] = S[t-1] + delta_S[t]
    E[t] = E[t-1] + delta_E[t]
    I[t] = I[t-1] + delta_I[t]
    R[t] = R[t-1] + delta_R[t]
    # observation process
    # assume no delays: a fraction of infectious exposures will be detected
    delta_severe[t] = delta_infective_exposures[t] * prop_severe;             # eq6
    severe_mean[t] = delta_severe[t]*pop                                      # eq7
    severe_obs_log[t] = rnorm(n=1,mean = log( severe_mean[t] ),sd = sigma_m ) # eq8
    severe_obs[t] = exp(severe_obs_log[t])
  }
  
  tibble(t_v,severe_obs,severe_mean,severe_obs_log,pop,S,I,E,R,delta_S,delta_E,delta_I,delta_R,delta_infective_exposures) -> df_hosp
  if (F){
    df_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
      geom_point(alpha=0.2)+
      geom_line(aes(y=delta_severe*pop),col="blue")
  }
  return(df_hosp)
}

# sim_wave_multi_voc