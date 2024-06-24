model_SIR_simple = function( params=NULL, dat=NULL, country_short_input, date_v_fit ){
  
  scenario_tag = "A"
  
  # ---- |-Fitting ----
  # prepare data for stan fit
  data_mock = dat %>% 
    filter(country_short == country_short_input, 
           target == params$SIR_simple$target, 
           agegroup == params$SIR_simple$agegroup) # 
  data_mock %<>% filter( date%in%date_v_fit )
  data_mock %>% ggplot(aes(date,value)) + geom_line()
  data_mock_fit = data_mock
  # Projection dates
  data_mock_project = data_mock
  data_mock_project$date = data_mock_project$date+365
  
  #mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = as.integer(data_mock_fit$value),
    n_week_project = nrow(data_mock_project),
    pop = 9e6,
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious
  )
  fit02=rstan::stan(
    file='./stan/SIR_simple.stan',
    chains=8 ,thin=8,iter=400,
    seed=12, cores = getOption("mc.cores", 1L),
    control=list(
      #adapt_delta=0.9,
      #max_treedepth=14
    ),
    data=stan_list
  ) # X mins
  
  # ---- |-Simulations ----
  #  inputs
  #- for each axis a list with available IDs
  
  ## transmission
  df_out = fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>%
    ungroup() %>%
    left_join(data_mock %>% select(date,country_short,agegroup,target) %>% mutate(t_vw = 1:n()), by="t_vw") %>%
    mutate(model = "SIR_simple",
           country_short = country_short,
           agegroup = agegroup,
           target = target,
           scenario_tag = scenario_tag,
           prediction_type = "sample"
    ) %>%
    rename(value = .value,
           sample_or_quantile = .draw) %>%
    select(-.chain, -.iteration, -.variable, -t_vw) 
  df_out %>% select(sample_or_quantile,date,value) %>% group_by(sample_or_quantile) %>% 
    nest() %>% rename(id=sample_or_quantile)  -> df
  list_transmission = df_to_list(df)
  list_transmission %>% names()
  ## vaccine uptake
  list_vaccine_id = list()
  list_vaccine_id[[1]] = list(
    vaccine_uptake = matrix(0*list_transmission[[1]]$value,ncol=1),
    VE_severe = 0.7
  )
  ## severity
  list_severe_id = list()
  list_severe_id[[1]] = list(
    
  )
  ## select needed IDs within each list & create mapping of the IDs across axes -> axis_ids_simulate
  axis_ids_simulate = tibble( round_id="2024-01",
                              transmission_id= c(1:100) ,
                              vaccine_id = 1,
                              severe_id = 1) %>% 
    ungroup() %>% mutate(sim_sample_id=1:n())
  axis_ids_simulate$sim = list(0)
  
  ## loop through axis_ids_simulate 
  for (sim_i in 1:nrow(axis_ids_simulate) ) {
    ## ---- |-Get axis IDs ----
    transmission_id = axis_ids_simulate$transmission_id[sim_i]
    vaccine_id = axis_ids_simulate$vaccine_id[sim_i]
    severe_id = axis_ids_simulate$severe_id[sim_i]
    
    ## ---- |-Prepare the details for each axis ----
    # prepare transmission
    transmission_df = list_transmission[[ transmission_id ]]
    incident_infections = transmission_df %>% select(value) %>% as.matrix() # format: [t,a]
    # prepare vaccination ( using list_vaccine_id )
    vaccine_uptake = list_vaccine_id[[vaccine_id]]$vaccine_uptake
    VE_severe = list_vaccine_id[[vaccine_id]]$VE_severe
    vax_sev = vaccine_severity_nowane(vaccine_uptake,VE_severe)
    # prepare natural severity ( using list_severe_id )
    severity_options = list_severe_id[[severe_id]]
    nat_sev = natural_severity(incident_infections,severity_options)
    severity_baseline = 0.5
    
    # run: transmission-to-severity 
    sev_fact = severity_factor(
      incident_infections,
      severity_baseline,
      severity_factor_vaccines=vax_sev$severity_factor_vaccines,
      severity_factor_natural=nat_sev$severity_factor_natural,
      severity_options
    )
    # run: combine all targets -> mysim
    mysim = combine_all_targets_SIR_simple(date_v=transmission_df$date,
                                           incident_infections,
                                           vaccine_uptake,
                                           incident_severe=sev_fact$incident_severe)
    axis_ids_simulate$sim[sim_i] = nest(mysim)[[1,1]]
    
    
    # Columns in the resulting df_out (as per: https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit):
    # model_id
    # round_id ["2024_2025_1_FLU1"]
    # scenario_id ["A","B"], target [allowed targets], location ["DE","FR"] 
    # pop_group ["0-12","13-65"], horizon [week integer], target_end_date [Date string ('YYYY-MM-DD')]
    # output_type ["sample"], output_type_id [string: "1","2","3",...], value [float limited to 2 decimals]
    
    # making up a fictive result as placeholder
    mdf = tibble(
      model_id = "ECDC_lefluflu",
      round_id = params$scenario_round_id,
      scenario_id = scenario_tag,
      target = "inc infection",
      location = country_short_input,
      pop_group = params$SIR_simple$agegroup,
      horizon = c(1,2,3,4),
      target_end_date = today()+(c(1,2,3,4)-1)*7,
      output_type = "sample",
      output_type_id = 1,
      value = c(12,25,12,12)
    )
    
    return(mdf)
  }
  
  
  if (F){ # support debugging
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe,
                           pop_infect
    ) %>% 
      mean_qi() -> xp; xp
    
    axis_ids_simulate %>% unnest(cols=sim) %>% 
      ggplot(aes(date,inc_death,group=sim_sample_id)) + geom_line()
    
    round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
    (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
    (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
    (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>% 
      mean_qi() %>% left_join(data_mock_fit %>% mutate(t_vw = 1:n()),by="t_vw") %>%
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      geom_point(aes(y=value),col="black") +
      labs(subtitle = paste("Austria 2022/2023: fit |",
                            "prob_severe:", mprob_severe,"\n",
                            "| S_ini:",mS_ini,
                            "| I_ini:",mI_ini,
                            "| prop inf:",mProp_inf)) -> p_cf0; p_cf0
    
    
  }
  
  return(axis_ids_simulate)
}

model_SIR_multiseason = function( params=NULL, 
                                  all_season=NULL, 
                                  target_input=NULL, 
                                  country_short_input , 
                                  pop_country ){
  
  # ---- |-Filtering and computing compound indicators ----
  all_season %>% 
    filter(country_short == country_short_input) %>% 
    select(-typing_sentinel,-typing_nonsentinel,-typing_combined) %>% 
    unnest(inc_iliari) %>% 
    filter(target==params$SIR_simple$target) -> all_season_filtered
  # filter:age groups
  all_season_filtered %>% 
    filter(agegroup==params$SIR_simple$agegroup) -> all_season_filtered
  # filter:pandemic seasons
  all_season_filtered %>% 
    filter(!season%in%params$SIR_multiseason$seasons_exclude) %>% 
    select(country_short,season,date,value) %>% 
    mutate(n=1:n()) -> all_season_fit
  
  if (target_input=="ili_typing_sentinel") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_sentinel)) %>% filter(indicator=="positivity") %>% 
      filter(!season%in%params$SIR_multiseason$seasons_exclude) %>% 
      select(date,value_typing=value)
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing/100)) %>% select(-value_typing)
  }
  if (target_input=="ili_typing_all") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter(!season%in%params$SIR_multiseason$seasons_exclude) %>% 
      select(date,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # ---- |-Project df and weekly to daily ----
  # project df
  start_year =year(today())
  season     = paste0(start_year,"/",start_year+1)
  start_date = ymd(paste0(start_year,params$season_start_monthday))
  end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
  date_v = seq(from=start_date,to=end_date,by="day")
  date_v_wed = date_v[weekdays(date_v)=="Wednesday"]
  all_season_project = tibble(country_short=country_short_input,
                              season=season,
                              date=date_v_wed,
                              value=NA)
  # from weekly, make daily
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_fit$date,season=all_season_fit$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) %>% 
    group_by(season) %>%  mutate(h=1:n(),
                                 season_start=case_when(h==1~1,h==2~2,.default=0) ) %>% 
    ungroup() %>% select(-h)-> all_season_fit_daily
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_project$date,season=all_season_project$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) -> all_season_project_daily
  # plotting
  p1 = all_season_fit %>% ggplot(aes(date,value)) + geom_line() + geom_rug()
  
  df_scenarios = 
    tibble(
      n=1:6,
      scenario_id=c("A","B","C","D","E","F"),
      axis_transmission=c(),
      axis_vax=c(1,1,2,2,3,3),
      axis_transmission=c(1,2,1,2,1,2)
    ) %>% left_join(
      tibble(axis_vax=c(1,2,3),axis_vax_name=c("opti","pess","null")) , by = join_by(axis_vax)
    ) %>% left_join(
      tibble(axis_transmission=c(1,2),axis_transmission_name=c("opti","pess")), by = join_by(axis_transmission)
    )
  
  # ---- |-Stan list and fit----
  stan_list = list(
    ## EXTRA stuff good to carry forward
    all_season_fit=all_season_fit,
    all_season_project=all_season_project,
    season_id_raw = fct_inorder(all_season_fit_daily$season) %>% levels() %>% enframe(),
    df_scenarios = df_scenarios,
    ## data relevated for the fit
    n_season = n_distinct(all_season_fit$season),
    n_week_fit = nrow(all_season_fit),
    n_day_fit = nrow(all_season_fit_daily),
    #
    n_age_groups = 2,
    #
    severe_obs_fit = ( all_season_fit %>% select(value) %>% mutate(value=value/2) %>% 
                         mutate(value=replace_na(value,0) %>% as.integer(),val=value) ),
    severe_obs_notna = as.integer(!is.na(all_season_fit$value)),
    season_start = as.integer(all_season_fit_daily$season_start),
    season_id = fct_inorder(all_season_fit_daily$season) %>% as.integer(),
    #
    pop = pop_country,
    #
    pop_age_group=matrix(c(pop_country/2,pop_country/2) ,nrow=2,ncol=1),
    contact_matrix=matrix(data=c(1,1,1,1),nrow=2,ncol=2),
    delta_vax=tibble( val1=rep(0,nrow(all_season_fit_daily)),
                      val2=rep(0,nrow(all_season_fit_daily))),
    # data relevant for projected scenarios
    n_week_project = nrow(all_season_project),
    n_day_project= nrow(all_season_project)*7,
    n_scenario = nrow(df_scenarios),
    axis_transmission = df_scenarios$axis_transmission,
    axis_vax = df_scenarios$axis_vax,
    delta_vax_opti=tibble( val1=rep(0,nrow(all_season_project_daily)),val2=rep(0,nrow(all_season_project_daily)) ),
    delta_vax_pess=tibble( val1=rep(0,nrow(all_season_project_daily)),val2=rep(0,nrow(all_season_project_daily)) ),
    delta_vax_null=tibble( val1=rep(0,nrow(all_season_project_daily)),val2=rep(0,nrow(all_season_project_daily)) ),
    # epi parameters
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious,
    ve_inf = params$ve_inf,
    ve_susc = params$ve_susc,
    ve_severe = params$ve_severe
  )
  ### make it 1 age group
  if (params$debug==TRUE) {
    stan_list$n_age_groups = 1
    stan_list$contact_matrix = matrix(data=c(1),nrow=1,ncol=1)
    stan_list$pop_age_group = matrix(c(pop_country) ,nrow=1,ncol=1)
    stan_list$severe_obs_fit =  ( all_season_fit %>% select(value) %>% 
                                    mutate(value=replace_na(value,0) %>% as.integer()) )
    stan_list$delta_vax = tibble( val1=rep(0,nrow(all_season_fit_daily)) )
    stan_list$delta_vax_opti = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_pess = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_null = tibble( val1=rep(0,nrow(all_season_project_daily)) )
  }
  
  p2 = NULL
  path_fit = paste0("../Big data/multiseason_age_vax",target_input,country_short_input,".Rdata")
  path_fit = "../Big data/multiseason_age_vaxili_typing_sentinelAT.Rdata"
  pr=paste("> Fitting:",target_input,"for",country_short_input,"... "); cleancat(green(pr))
  if (params$debug==F) {
    
    fit00=rstan::stan(
      file='./stan/SIR_multiseason_age_vax.stan',
      chains=1 ,thin=1,iter=100,
      seed=12, cores = getOption("mc.cores", 1L),
      control=list(
        adapt_delta=0.95
        #max_treedepth=14
      ),
      data=stan_list
    ) # 2.3 mins
   
    
    save(fit00,stan_list,file = path_fit)
  } else {
    load(path_fit)
  }
  pr=paste("> Fitting:",target_input,"for",country_short_input,"Done \n"); cleancat(green(pr))
  
  # ---- |-Extract parameters and plot ----
  if (F){
    precis(fit00,pars=c("SIR_ini_mu"),depth = 3)
    precis(fit00,pars=c("Rnull_eff"),depth = 2)
    precis(fit00,pars=c("prop_severe"),depth = 2)
    precis(fit00,pars=c("sigma_i"),depth = 2) # 3.66
    precis(fit00,pars=c("SIR_ini"),depth = 3)
  }
  
  
  p2 = fit00 %>% gather_draws(delta_severe_abs_weekly_sum[n]) %>% 
    filter(.draw%in%c(1:20)) %>% select(-.chain,-.iteration) %>% ungroup() %>% 
    right_join(all_season_fit,by = join_by(n)) %>% 
    ggplot(aes(date,value)) + 
    geom_line(aes(y=.value,group=.draw),col="lightblue") + geom_line() + coord_cartesian(ylim=c(0,10000))
  
  # ---- |-Compile output ----
  modl = 
    list(
      fit = fit00,
      stan_list = stan_list,
      pdata = p1,
      pfit = p2
    )
  
  return(modl)
}

model_SIR_simple_r0 = function( params=NULL, all_season=NULL, target_input=NULL,pop_country, country_short_input, date_v_fit,season){
  
  # ---- |-Filtering ----
  all_season %>% 
    filter(country_short == country_short_input) %>% 
    select(-typing_sentinel,-typing_nonsentinel,-typing_combined) %>% 
    unnest(inc_iliari) %>% 
    filter(target==params$SIR_simple$target) -> all_season_filtered
  # filter:age groups
  all_season_filtered %>% 
    filter(agegroup==params$SIR_simple$agegroup) -> all_season_filtered
  # filter:time
  all_season_filtered %>% 
    filter( date%in%date_v_fit ) %>% 
    select(country_short,season,date,value) %>% 
    mutate(n=1:n()) -> all_season_fit
  
  if (target_input=="ili_typing_sentinel") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_sentinel)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value)
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing/100)) %>% select(-value_typing)
  }
  if (target_input=="ili_typing_all") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # prepare data for stan fit
  all_season_fit %>% ggplot(aes(date,value)) + geom_line()
  data_mock_fit = all_season_fit
  # Projection dates
  data_mock_project = all_season_fit
  data_mock_project$date = data_mock_project$date+365
  
  #mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = as.integer(replace_na(data_mock_fit$value,0)),
    severe_obs_notna = as.integer(!is.na(data_mock_fit$value)),
    n_week_project = nrow(data_mock_project),
    pop = pop_country,
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious
  )
  fit02=rstan::stan(
    file='./stan/SIR_simple_nas.stan',
    #chains=1 ,thin=1,iter=300,
    chains=6 ,thin=6,iter=1500,
    seed=12, cores = getOption("mc.cores", 1L),
    control=list(
      #adapt_delta=0.9,
      #max_treedepth=14
    ),
    data=stan_list
  ) # X mins
  est_Rnull = precis(fit02,pars = c("Rnull_eff"))
  
  mout = tibble(
    country_short = country_short_input,
    season = season,
    Rnull = est_Rnull$result[1],
    Rnull_Rhat = est_Rnull$result[6]
  )
  
  
  if (F){ # plotting to support debugging
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe,
                           pop_infect
    ) %>% 
      mean_qi() -> xp; xp
    
    
    round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
    (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
    (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
    (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>% 
      mean_qi() %>% left_join(data_mock_fit %>% mutate(t_vw = 1:n()),by="t_vw") %>%
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      geom_point(aes(y=value),col="black") +
      labs(subtitle = paste("Austria 2022/2023: fit |",
                            "prob_severe:", mprob_severe,"\n",
                            "| S_ini:",mS_ini,
                            "| I_ini:",mI_ini,
                            "| prop inf:",mProp_inf)) -> p_cf0; p_cf0
    
    
  }
  return(mout)
}