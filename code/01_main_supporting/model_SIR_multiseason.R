model_SIR_multiseason = function( params=NULL, 
                                  all_season=NULL, 
                                  target_input=NULL, 
                                  country_short_input , 
                                  pop_country, 
                                  vax_country,
                                  contacts){
  
  # ---- |-Obtain the fitting dataframe from the data ----
  # with age groups in separate columns
  all_season_fit_wide = wrangle_fit_df(all_season,params,country_short_input,target_input)
  
  # ---- |-Project df ----
  start_year = year(today())
  season     = paste0(start_year,"/",start_year+1)
  start_date = ymd(paste0(start_year,params$season_start_monthday))
  end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
  date_v = seq(from=start_date,to=end_date,by="day")
  date_v_wed = date_v[weekdays(date_v)=="Wednesday"][1:52]
  date_v_mon = date_v[weekdays(date_v)=="Monday"][1:52]
  all_season_project = tibble(country_short=country_short_input,
                              season=season,
                              date_mon=date_v_mon,
                              date_wed=date_v_wed,
                              value=NA) %>% 
    mutate(week_id=1:n()) %>% 
    left_join(data$helpers_respicompass$iso_weeks,by=c("date_mon"="start_week_day"))
  
  data$helpers_respicompass$iso_weeks
  # make daily version of the data frame - from some daily indicators that the model needs
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_fit_wide$date,season=all_season_fit_wide$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) %>% 
    group_by(season) %>%  mutate(h=1:n(),
                                 season_start=case_when(h==1~1,h==2~2,.default=0) ) %>% 
    ungroup() %>% select(-h) -> all_season_fit_daily
  
  crossing( country_short=country_short_input,
            nesting(data_w=all_season_project$date_mon,season=all_season_project$season), 
            d_shift=c(-6:0) ) %>% 
    mutate(date=data_w+d_shift) %>% select(country_short,season,date) -> all_season_project_daily
  # plotting
  p1 = all_season_fit_wide %>% ggplot(aes(date,age_total)) + geom_line() + geom_rug()
  
  
  # helpers for stan list
  df_scenarios = params$scenarios
  
  pop_pyramid = data$demography_respicast$population_pyramid %>% filter(country==EU_long(country_short_input))
  pop_age_group = pop_pyramid$population 
  age_groups = params$SIR_multiseason$age_groups
  n_age_groups = length(age_groups)
  df_agegroups_ecdc = tibble(agegroup_id=c(0:n_age_groups),
                             age_group_ecdc=c("age_total", age_groups) )
  df_age_translate = tibble(age_group_ecdc=df_agegroups_ecdc$age_group_ecdc,
                            age_group_respicompass=c("total","0-4","5-14","15-64","65+"))
  df_agegroups = df_agegroups_ecdc %>% left_join(df_age_translate,by = join_by(age_group_ecdc))
  z_proj = rep(0,nrow(all_season_project_daily))
  z_fit  = rep(0,nrow(all_season_fit_daily))
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
    n_age_groups = n_age_groups,
    #
    ili_obs_fit = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ replace_na(.,0) ) %>% mutate_all(~as.integer(.) ),
    ili_obs_notna = all_season_fit_wide %>% 
      select( any_of(params$SIR_multiseason$age_groups) ) %>% 
      mutate_all(~ !is.na(.) ) %>% mutate_all(~as.integer(.) ),
    #
    season_start = as.integer(all_season_fit_daily$season_start),
    season_id = fct_inorder(all_season_fit_daily$season) %>% as.integer(),
    #
    pop = sum(pop_age_group), 
    #
    pop_age_group=matrix(pop_age_group ,nrow=n_age_groups,ncol=1),
    #contact_matrix=matrix(data=rep(1/n_age_groups,n_age_groups^2),nrow=n_age_groups,ncol=n_age_groups),
    contact_matrix=matrix(data=rep(1,n_age_groups^2),nrow=n_age_groups,ncol=n_age_groups),
    delta_vax=tibble( A=z_fit,B=z_fit,C=z_fit,D=z_fit) %>% mnaming(age_groups),
    # data relevant for projected scenarios
    n_week_project = nrow(all_season_project),
    n_day_project= nrow(all_season_project)*7,
    n_scenario = nrow(df_scenarios),
    axis_transmission = df_scenarios$axis_transmission,
    axis_vax = df_scenarios$axis_vax,
    delta_vax_real=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_opti=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_pess=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    delta_vax_null=tibble( A=z_proj,B=z_proj,C=z_proj,D=z_proj) %>% mnaming(age_groups),
    # epi parameters
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious,
    ve_spread = params$ve_spread,
    ve_inf = params$ve_inf,
    ve_ili_cond_inf = params$ve_ili_cond_inf
  )
  # Add vaccination to Oct 1st to the second age group
  ind_vax = which(date_v == paste0(year(min(date_v)),"-10-01"))
  stan_list$delta_vax_real$age_65_99[ind_vax] = (vax_country$higher_vax_coverage + vax_country$lower_vax_coverage)/2
  stan_list$delta_vax_opti$age_65_99[ind_vax] = vax_country$higher_vax_coverage
  stan_list$delta_vax_pess$age_65_99[ind_vax] = vax_country$lower_vax_coverage
  stan_list$delta_vax_null$age_65_99[ind_vax] = 0
  
  
  ############## Add artifically generated ili values ###############
  
  stan_list = generate_ili_epi_test(stan_list)
  
  
  ###################################################################
  
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
    stan_list$delta_vax_pess = tibble( val1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_null = tibble( val1=rep(0,nrow(all_season_project_daily)) )
  }
  
  p2 = NULL
  path_fit = paste0("../Big data/multiseason_age_vax",target_input,country_short_input,".Rdata")
  pr=paste("> Now fitting:",target_input,"for",country_short_input,"... "); cleancat(green(pr))
  if (T) {
    
    # using a preious fit, set initial conditions
    ini_tune = F
    if (ini_tune==T) {
      load(path_fit) # loading fit00
      fit_means = get_posterior_mean(fit00) # extract mean estimates
      row.names(fit_means)[1:32] # print parameter names
      # mp="sigma_i";mcmc_areas(fit00,mp);precis(fit00,depth=3,mp)
      #
      myl = vector(mode = "list", length = 32)
      myl[1:32] = fit_means[1:32]
      names(myl) = row.names(fit_means)[1:32]
      init_fun = function(...) myl
      save(init_fun,file=paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata") )
      save_as_general = T
      if (save_as_general) save(init_fun,file="output/ini_SIR__multiseason_age_vax.Rdata")
    }
    if (ini_tune==F) {
      # a specific country run will look for a country-specfic initial contition function file, if that does not exist, it loads a general one
      if (file.exists(paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata"))){
        load(file=paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata"))  
      } else{
        load(file="output/ini_SIR__multiseason_age_vax.Rdata")
      }
    }
    
    fit00=rstan::stan(
      file='./stan/SIR_multiseason_age_vax.stan',
      #chains=1 ,thin=1,iter=150, # a "debug run"
      chains=4, thin=4, iter=1000, # a "long run" 
      seed=13, cores = getOption("mc.cores", 1L),
      control=list(
        # adapt_delta=0.95, # look into increasing this, 0.98 or 0.99
        # max_treedepth=10 # look into increasing this to, 15, 20 ect
      ),
      data=stan_list
      #init = init_fun
    ) # 8.5 hrs, run time will scale with  inter, and is a function of adapt_delta and max_treedepth, and is a function of luck
    
    
    save(fit00,stan_list,file = path_fit)
  } else {
    load(path_fit)
  }
  pr=paste("> Fitting:",target_input,"for",country_short_input,"Done \n"); cleancat(green(pr))
  
  # ---- |-Block to sense check etc, Extract parameters and plot, see convergence ----
  if (F){
    # plot without Rethinking package
    fit00@model_pars # see which parameters are there
    # see parameter names with dimensions
    fit_means = get_posterior_mean(fit00) # extract mean estimates
    row.names(fit_means)[1:32]
    mp="SIR_ini[3,1,1]"; mcmc_areas(fit00,mp)
    mp="reciprocal_phi"; summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary
    
    
    # https://rstudio.github.io/cheatsheets/bayesplot.pdf
    # for quick convergence check: n_eff Rhat ( see other packages than Rethinking)
    # using Rethinkingp akcage
    
    # precis(fit00,pars=c("SIR_ini_mu"),depth = 3)
    # precis(fit00,pars=c("Rnull_eff"),depth = 2)
    # precis(fit00,pars=c("prop_severe"),depth = 2)
    # precis(fit00,pars=c("sigma_i"),depth = 2) # 3.66
    # precis(fit00,pars=c("SIR_ini"),depth = 3) 
  }
  
  # extract fit
  modelled_fit = fit00 %>% gather_draws(gen_ili_obs_fit_sum[n]) %>% 
    filter(.draw%in%c(1:20)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    right_join(all_season_fit_wide,by = join_by(n)) 
  modelled_fit %>% ggplot(aes(date,age_total)) + geom_line() + 
    geom_line(aes(y=.value,group=.draw),col="lightblue") +
    coord_cartesian(ylim = c(0,2*modelled_fit$age_total %>% max(na.rm=T)))
  
  # extract projections
  modelled_proj = fit00 %>% 
    gather_draws(gen_ili_u_obs_project[scen_id,week_id,agegroup_id], # if you want to apply changes here, do look up the useful gather_draws {tidybayes} syntax
                 gen_ili_v_obs_project[scen_id,week_id,agegroup_id],
                 gen_ili_u_obs_project_sum[scen_id,week_id],
                 gen_ili_v_obs_project_sum[scen_id,week_id]) %>% 
    filter(.draw%in%c(1:20)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% # remove unneeded columns and grouping
    mutate(vax_status=case_when(
      .variable=="gen_ili_u_obs_project"~"vaxNo",
      .variable=="gen_ili_v_obs_project"~"vaxYes",
      .variable=="gen_ili_u_obs_project_sum"~"vaxNo",
      .variable=="gen_ili_v_obs_project_sum"~"vaxYes"
    )) %>% 
    mutate(agegroup_id=case_when(
      .variable=="gen_ili_u_obs_project_sum"~0,
      .variable=="gen_ili_v_obs_project_sum"~0,
      TRUE ~ agegroup_id
    )) %>% 
    left_join(df_scenarios,by = join_by(scen_id)) %>% # add scenario info
    left_join(df_agegroups,by = join_by(agegroup_id)) %>% # add agegroup info
    left_join(all_season_project, by=join_by(week_id)) %>% # add week info
    mutate(model_id=params$scenario_model,round_id=params$scenario_round_id, # add needed columns
           target="ili_plus",output_type="sample") %>% 
    unite(col="pop_group",age_group_respicompass,vax_status,sep="_") %>% 
    select( round_id, # select according to submission definition
            scenario_id,
            target=target,
            location=country_short,
            pop_group=pop_group,
            horizon=horizon,
            target_end_date=end_week_day,
            date_mon=date_mon, # keeping this in to have dates for all weeks beyond data$helpers_respicompass$iso_weeks
            output_type=output_type,
            output_type_id=.draw,
            value=.value) 
  # Required columns in df_for_submission (as per: https://github.com/european-modelling-hubs/RespiCompass/wiki/Submission-format):
  
  # ---- |-Compile output ----
  modl = 
    list(
      stan_list = stan_list,
      modelled_fit = modelled_fit,
      modelled_proj = modelled_proj
    )
  
  return(modl)
}

# functions supporting model_SIR_multiseason()
wrangle_fit_df = function(all_season,params,country_short_input,target_input){
  
  if (target_input=="respicompass_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(respicompass_ili_plus) %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value)-> all_season_fit
  }
  
  if (target_input=="ili") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(inc_iliari) %>% 
      filter(target=="ILIconsultationrate") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value) -> all_season_fit
  }
  
  if (target_input=="own_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(inc_iliari) %>% 
      filter(target=="ILIconsultationrate") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value) -> all_season_ili
    
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # take care of age groups
  all_season_fit_wide = all_season_fit %>% 
    pivot_wider(names_from = agegroup, values_from = value) %>% 
    mutate(n=1:n())
  
  return(all_season_fit_wide)
} 
