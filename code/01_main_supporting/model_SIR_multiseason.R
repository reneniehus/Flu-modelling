fit_with_stan = function(params,stan_list,mod_path,all_season_fit_wide) {
  
  # initiate output list
  mout=list()
  
  # compile the stan model 
  m <- stan_model(file=mod_path)
  browser()
  # run the model fit
  if (T) {
    start_time <- Sys.time()
    fit00=rstan::vb(
      m,
      algorithm = "meanfield", # variational inference algorithm
      grad_samples=5 , # samples to determine the gradient ( 2 is slower than 5, )
      tol_rel_obj = 0.01,
      iter=1000, # 
      output_samples = 500,
      #chains=2, thin=2, iter=300, # a "long run" 
      seed=12, # seed for pseudo-random numbers to ensure reproducibility
      data=stan_list # data input into the model
    ) # 
    end_time <- Sys.time(); end_time - start_time
  } # 1.6 min
  
  browser()
  # plot the fit against fitted data
  modelled_fit = fit00 %>% gather_draws(gen_ili_obs_fit_sum[n]) %>% 
    filter(.draw%in%c(1:20)) %>% # filter a number of posterior draws
    select(-.chain,-.iteration) %>% ungroup() %>% 
    right_join( all_season_fit_wide %>% mutate(age_total = stan_list$ili_obs_fit$age_1 ),
                by = join_by(n)) %>%
    group_by(date) %>%
    mutate(mean_value = mean(.value)) %>%
    ungroup()
  modelled_proj = fit00 %>% gather_draws(gen_ili_t_obs_project_sum[scen,week_id]) %>% 
    filter(.draw%in%c(1:50)) %>%
    select(-.chain,-.iteration) %>% ungroup() %>% 
    right_join( stan_list$all_season_project ,
                by = join_by(week_id)) %>% 
    group_by(date=date_wed,scen) %>%
    mutate(mean_value = mean(.value)) %>%
    ungroup()
    
  p1 = modelled_fit %>% ggplot(aes(date,age_total)) + geom_line() + 
    geom_line(aes(y=mean_value),col="lightblue") +
    geom_line(data=modelled_proj,aes(col=as.factor(scen),y=mean_value)) +
    coord_cartesian(ylim = c(0,2*modelled_fit$age_total %>% max(na.rm=T))) ; p1
  
  library(rethinking)
  precis(fit00, pars="prop_ili_mu",depth=2) # prop_ili_mu
  stan_list$cum_ili_obs_log
  
  if (F) source("code/01_main_supporting/old_stan_fit_code.R")
  
  # add into output list 
  mout$stan_list = stan_list
  mout$fit = fit00
  mout$modelled_proj = extract_projections(params,fit00,n_iter=20,
                                           stan_list$df_scenarios,
                                           stan_list$df_agegroups,
                                           stan_list$all_season_project)
  mout$plot_fit = p1
  return(mout)
}

# fitting with sequential ABC 
fit_with_eabc = function(params,stan_list,mod_path) {
  
  # Define priors
  myPriors <- list('S01' = c("unif",0,0.5),
                   'S02' = c("unif",0,0.5),
                   'S03' = c("unif",0,0.5))
  myPriors <- list('S01' = c("unif",0,0.5))
  
  if (T){
    tic()
    x<-generate_ili_epi_test( c(1,0.3),stan_list) 
    toc()
  }
  
  # Wrap up model in function that outputs summary stats
  myModel <- function(par){
    
    stan_list_f = generate_ili_epi_test(par,stan_list)
    
    return( stan_list_f$ili_obs_fit$age_1[stan_list_f$ili_obs_notna$age_1==1] ) 
  }
  
  # Define targets 
  myTarget <- c( stan_list$ili_obs_fit$age_1[stan_list$ili_obs_notna$age_1==1] )
  
  # 
  dist_euc <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
  dist_euc(myTarget,myModel( c(1,0.1)))
  
  # Run ABC-SMC (this should be parallelised, see package help)
  library(EasyABC)
  
  rval <- ABC_sequential(method = "Beaumont", 
                         model = myModel, 
                         prior = myPriors, 
                         nb_simul = 5, 
                         summary_stat_target = myTarget,
                         #n_cluster=8,
                         tolerance_tab = c(188521,188520),
                         use_seed=TRUE,
                         progress_bar=T
  )
  
  # Plot posteriors
  hist(rval$param[,1])
  hist(rval$param[,2])
  plot(rval$param[,1], rval$param[,2])
  
}

# computing the data frame from all_season_country
wrangle_fit_df = function(params,data,all_season_country,country_short_input,target_input){
  
  all_season = all_season_country
  
  if (target_input=="respicompass_ili_plus") {
    all_season %>% 
      filter(country_short == country_short_input) %>% 
      unnest(respicompass_ili_plus) %>% 
      filter(season%in%params$SIR_multiseason$seasons_include) %>% 
      select(country_short,date,season,agegroup,value)-> all_season_fit
    # impute summer low-activity
    all_season_fit %>% mutate(
      value=ifelse( date%in%params$summer_low_dates, 0 , value )
    ) -> all_season_fit
    
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

# computing a list with all input required by the model
make_stan_list = function(params,data,all_season_fit_wide,country_short_input,vax_country,pop_country){
  # helpers for the dataframes
  start_year = year(today())
  season     = paste0(start_year,"/",start_year+1)
  start_date = ymd(paste0(start_year,params$season_start_monthday))
  end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
  date_v = seq(from=start_date,to=end_date,by="day")
  date_v_wed = date_v[weekdays(date_v)=="Wednesday"][1:52]
  date_v_mon = date_v[weekdays(date_v)=="Monday"][1:52]
  
  # dataframe for projections
  all_season_project = tibble(country_short=country_short_input,
                              season=season,
                              date_mon=date_v_mon,
                              date_wed=date_v_wed,
                              value=NA) %>% 
    mutate(week_id=1:n()) %>% 
    left_join(data$helpers_respicompass$iso_weeks,by=c("date_mon"="start_week_day"))
  
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
    df_agegroups = df_agegroups,
    ili_obs_fit_date = all_season_fit_wide$date,
    ili_obs_project = all_season_project$date_wed,
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
    season_id_day = fct_inorder(all_season_fit_daily$season) %>% as.integer(),
    season_id_week = fct_inorder(all_season_fit_wide$season) %>% as.integer(),
    #
    pop = sum(pop_age_group), 
    #
    pop_age_group=matrix(pop_age_group ,nrow=n_age_groups,ncol=1),
    #contact_matrix=matrix(data=rep(1/n_age_groups,n_age_groups^2),nrow=n_age_groups,ncol=n_age_groups),
    contact_matrix=matrix(data=rep(1/n_age_groups,n_age_groups^2),nrow=n_age_groups,ncol=n_age_groups),
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
    ve_ili_cond_inf = params$ve_ili_cond_inf,
    # daily steps
    n_daily_time_steps = 1,
    # priors
    sigma_cum_ili = 5,
    prior_sigma_prop_ili = 2,
    prior_sigma_i = 5,
    prior_sigma_s = 1
  )
  # Add vaccination to Oct 1st to the second age group
  ind_vax = which(date_v == paste0(year(min(date_v)),"-10-01"))
  stan_list$delta_vax_real$age_65_99[ind_vax] = (vax_country$higher_vax_coverage + vax_country$lower_vax_coverage)/2
  stan_list$delta_vax_opti$age_65_99[ind_vax] = vax_country$higher_vax_coverage
  stan_list$delta_vax_pess$age_65_99[ind_vax] = vax_country$lower_vax_coverage
  stan_list$delta_vax_null$age_65_99[ind_vax] = 0
  stan_list$delta_vax$age_65_99[ind_vax] = (vax_country$higher_vax_coverage + vax_country$lower_vax_coverage)/2
  stan_list$daily_counter_fit = rep(1:stan_list$n_day_fit, each=stan_list$n_daily_time_steps)
  stan_list$daily_counter_proj = rep(1:stan_list$n_day_project, each=stan_list$n_daily_time_steps)
  
  stan_list$daily_daystart_fit = rep(1:stan_list$n_daily_time_steps, each=stan_list$n_day_fit)
  stan_list$daily_daystart_proj = rep(1:stan_list$n_daily_time_steps, each=stan_list$n_day_proj)
  
  # summary targets
  stan_list$cum_ili_obs_log = rowsum(x=stan_list$ili_obs_fit,group=stan_list$season_id_week,na.rm = T) %>% rowSums() %>% log()
  stan_list$n_ili_obs_notna = rowsum(x=stan_list$ili_obs_notna,group=stan_list$season_id_week,na.rm = T) %>% rowSums()
  stan_list$weight_obs_epi =  1.0 #1/mean( stan_list$n_ili_obs_notna ) 
  stan_list$weight_cum_ili =  1.0 # 
  
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
      age_1=as.integer((age_00_04+age_05_14)==2),
      age_2=as.integer((age_15_64+age_65_99)==2)
    )
    stan_list$delta_vax_real = tibble( age_1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_opti = tibble( age_1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_pess = tibble( age_1=rep(0,nrow(all_season_project_daily)) )
    stan_list$delta_vax_null = tibble( age_1=rep(0,nrow(all_season_project_daily)) )
  }
  ### for debugging: make it 1 age group
  if (T) {
    stan_list$n_age_groups = 1
    stan_list$contact_matrix = matrix(data=c(1),nrow=1,ncol=1)
    stan_list$pop_age_group = matrix(c(pop_country) ,nrow=1,ncol=1)
    stan_list$ili_obs_notna =  stan_list$ili_obs_notna %>% transmute(
      age_1=as.integer( (age_00_04+age_05_14+age_15_64+age_65_99) == 4)
    )
    stan_list$ili_obs_fit =  stan_list$ili_obs_fit %>% transmute(
      age_1= age_00_04+age_05_14+age_15_64+age_65_99
    )
    stan_list$delta_vax_real = stan_list$delta_vax_real %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_opti = stan_list$delta_vax_opti %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_pess = stan_list$delta_vax_pess %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax_null = stan_list$delta_vax_null %>% select(age_65_99) %>% rename(age_1=age_65_99)
    stan_list$delta_vax = tibble( age_1=rep(0,nrow(all_season_fit_daily)) )
  }
  return(stan_list)
}

# extract projections from the model in the format required by RespiCompass
extract_projections = function(params,fit00,n_iter,df_scenarios,df_agegroups,all_season_project){
  # extract projections
  modelled_proj = fit00 %>% 
    gather_draws(gen_ili_u_obs_project[scen_id,week_id,agegroup_id], # if you want to apply changes here, do look up the useful gather_draws {tidybayes} syntax
                 gen_ili_v_obs_project[scen_id,week_id,agegroup_id],
                 gen_ili_u_obs_project_sum[scen_id,week_id],
                 gen_ili_v_obs_project_sum[scen_id,week_id]) %>% 
    filter(.draw%in%c(1:n_iter)) %>% # filter a number of posterior draws
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
  return(modelled_proj)
}
