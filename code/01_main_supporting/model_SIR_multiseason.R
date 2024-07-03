model_SIR_multiseason = function( params=NULL, 
                                  all_season=NULL, 
                                  target_input=NULL, 
                                  country_short_input , 
                                  pop_country, 
                                  vax_country){
  
  # ---- |-Filtering and computing compound indicators ----
  if (target_input=="respicompass_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      select(-typing_sentinel,-typing_nonsentinel,-typing_combined) %>% 
      unnest(respicompass_ili_plus) %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value)-> all_season_fit
  } else {
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
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
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
  }
  
  #
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
  
  
  # take care of age groups
  all_season_fit_wide = all_season_fit %>% 
    pivot_wider(names_from = agegroup, values_from = value) %>% 
    mutate(n=1:n())
  
  # make daily version of the data frame - from some daily indicators that the model needs
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_fit$date,season=all_season_fit$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) %>% 
    group_by(season) %>%  mutate(h=1:n(),
                                 season_start=case_when(h==1~1,h==2~2,.default=0) ) %>% 
    ungroup() %>% select(-h) -> all_season_fit_daily
  
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_project$date,season=all_season_project$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) -> all_season_project_daily
  # plotting
  p1 = all_season_fit %>% ggplot(aes(date,value,color=agegroup)) + geom_line() + geom_rug()
  
  # scenarios as per user settings
  df_scenarios = params$scenarios
  
  # ---- |-Stan list and fit----
  stan_list = list(
    ## EXTRA stuff good to carry forward
    all_season_fit=all_season_fit_wide,
    all_season_project=all_season_project,
    season_id_raw = fct_inorder(all_season_fit_daily$season) %>% levels() %>% enframe(),
    df_scenarios = df_scenarios,
    ## data relevated for the fit
    n_season = n_distinct(all_season_fit_wide$season),
    n_week_fit = nrow(all_season_fit_wide),
    n_day_fit = nrow(all_season_fit_daily),
    n_age_groups = length(params$SIR_multiseason$age_groups),
    #
    ili_obs_fit = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ replace_na(.,0) ) %>% mutate_all(~as.integer(.) ),
    ili_obs_notna = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ !is.na(.) ) %>% mutate_all(~as.integer(.) ),
    season_start = as.integer(all_season_fit_daily$season_start),
    season_id = fct_inorder(all_season_fit_daily$season) %>% as.integer(),
    #
    pop = sum(1315044,7789728), # hard coded for AT
    #
    pop_age_group=matrix(c(1315044,7789728) ,nrow=2,ncol=1),
    contact_matrix=matrix(data=c(1,1,1,1),nrow=2,ncol=2),
    delta_vax=tibble( val1=rep(0,nrow(all_season_fit_daily)),
                      val2=rep(0,nrow(all_season_fit_daily))),
    # data relevant for projected scenarios
    n_week_project = nrow(all_season_project),
    n_day_project= nrow(all_season_project)*7,
    n_scenario = nrow(df_scenarios),
    axis_transmission = df_scenarios$axis_transmission,
    axis_vax = df_scenarios$axis_vax,
    delta_vax_real=tibble( val1=rep(0, nrow(all_season_project_daily)), val2=rep(0,nrow(all_season_project_daily)) ),
    delta_vax_opti=tibble( val1=rep(0, nrow(all_season_project_daily)), val2=rep(0,nrow(all_season_project_daily)) ),
    delta_vax_pess=tibble( val1=rep(0, nrow(all_season_project_daily)), val2=rep(0,nrow(all_season_project_daily)) ),
    delta_vax_null=tibble( val1=rep(0, nrow(all_season_project_daily)), val2=rep(0,nrow(all_season_project_daily)) ),
    # epi parameters
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious,
    ve_spread = params$ve_spread,
    ve_inf = params$ve_inf,
    ve_ili_cond_inf = params$ve_ili_cond_inf
  )
  
  # Add vaccination to Oct 1st to the second age group
  ind_vax = which(date_v == paste0(year(min(date_v)),"-10-01"))
  stan_list$delta_vax_opti$val2[ind_vax] = vax_country$higher_vax_coverage
  stan_list$delta_vax_pess$val2[ind_vax] = vax_country$lower_vax_coverage
  stan_list$delta_vax_null$val2[ind_vax] = 0
  
  ### for debugging: make it 2 age groups
  if (F) {
    stan_list$n_age_groups = 2
    stan_list$contact_matrix = matrix(data=c(1,1,1,1),nrow=2,ncol=2)
    stan_list$pop_age_group = matrix(c(1315044,7789728) ,nrow=2,ncol=1)
    stan_list$ili_obs_fit =  all_season_fit_wide %>% transmute(
      age_1=replace_na(age_00_04+age_05_14,0) %>% as.integer(),
      age_2=replace_na(age_15_64+age_65_99,0) %>% as.integer()
    )
    stan_list$ili_obs_notna = all_season_fit_wide %>% transmute(
      age_1=as.integer(!is.na(age_00_04+age_05_14)),
      age_2=as.integer(!is.na(age_15_64+age_65_99))
    )
    stan_list$delta_vax_real = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_opti = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_opti = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_pess = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_null = tibble( val1=rep(0,nrow(all_season_project_daily)) )
  }
  ### for debugging: make it 1 age group
  if (F) {
    stan_list$n_age_groups = 1
    stan_list$contact_matrix = matrix(data=c(1),nrow=1,ncol=1)
    stan_list$pop_age_group = matrix(c(pop_country) ,nrow=1,ncol=1)
    stan_list$ili_obs_notna =  all_season_fit_wide %>% transmute(
      age_1=as.integer(!is.na(age_00_04+age_05_14+age_15_64+age_65_99))
    )
    stan_list$ili_obs_fit =  all_season_fit_wide %>% transmute(
      age_1=replace_na(age_00_04+age_05_14+age_15_64+age_65_99,0) %>% as.integer()
    )
    stan_list$delta_vax_real = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_opti = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_opti = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_pess = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_null = tibble( val1=rep(0,nrow(all_season_project_daily)) )
  }
  
  p2 = NULL
  path_fit = paste0("../Big data/multiseason_age_vax",target_input,country_short_input,".Rdata")
  path_fit = "../Big data/multiseason_age_vaxili_typing_sentinelAT.Rdata"
  pr=paste("> Fitting:",target_input,"for",country_short_input,"... "); cleancat(green(pr))
  if (T) {
    
    # using a preious fit, set initial conditions
    ini_tune = F
    if (ini_tune==T) {
      load(path_fit) # loading fit00
      fit_means = get_posterior_mean(fit00)
      mp="SIR_ini_mu[1,3]";mcmc_areas(fit00,mp);precis(fit00,depth=3,mp)
      #
      myl = vector(mode = "list", length = 24)
      myl[1:24] = fit_means[1:24]
      names(myl) = row.names(fit_means)[1:24]
      init_fun = function(...) myl
      save(init_fun,file="output/ini_SIR__multiseason_age_vax.Rdata")
    }
    if (ini_tune==F) {
      load(file="output/ini_SIR__multiseason_age_vax.Rdata")
    }
    
    fit00=rstan::stan(
      file='./stan/SIR_multiseason_age_vax.stan',
      chains=1 ,thin=1,iter=150,
      seed=12, cores = getOption("mc.cores", 1L),
      control=list(
        # adapt_delta=0.97,
        # max_treedepth=14
      ),
      data=stan_list,
      init = init_fun
    ) # X mins
    
    
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
  
  
  p2 = fit00 %>% gather_draws(delta_ili_abs_weekly_sum[n]) %>% 
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