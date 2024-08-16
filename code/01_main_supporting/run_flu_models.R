run_flu_models = function( params=NULL , data=NULL ){
  
  if (F){
    source("./code/00_main.R")
  }
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Initiating desired output list ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_out = list(
    time_of_execution = now(),    # time-stamp
    figs_prefit = NULL,           # figures of data prior to entering the fitting functions
    submission = NULL,     # for each model a clean dataframe following submission format (see above)
    output = NULL           # for each model additional output 
  )
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Running selected models ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if ( "SIR_simple_multi_season" %in% params$models_to_run ){
    start_time <- Sys.time()
    
    pr=paste("Initiating SIR_simple_multi_season \n"); cat(green(pr))
    
    # ---- |-Data for all countries ----
    all_season = data_into_all_season(data,params,withforce=F); df_out$figs_prefit$fit_seasons_countries <- plot_add_season(all_season)
    contacts = transform_contracts(data,params) # transform the contact matrixes for model requirements
    target_input_v = params$SIR_simple_multi_season$target
    country_short_input_v = all_season %>% filter_log(ili_plus_sum>0) %>% pull(country_short) %>% unique()
    
    modl = list()
    
    # ---- |-Run model for each country ----
    target_input=target_input_v[1]
    country_short_input_v = country_short_input_v # params$run_countries # c("IT","AT")
    for (country_short_input in country_short_input_v ) { # country_short_input="IT"
      # ---- |-Prepare country specific data ----
      pop_country = data$demography_respicast$population_pyramid %>% 
        filter(country==EU_long(country_short_input)) %>% pull(population) %>% sum()
      vax_country = data$vax$data_vax %>% filter( location_name == EU_long(country_short_input) ); if (nrow(vax_country) != 1) stop("Vaccination data is wrong format: either no data or too many rows")
      all_season_country = all_season %>% filter( country_short == country_short_input )
      pr=paste(target_input,"for",country_short_input,"with population:",pop_country,"\n"); cat(yellow(pr))
      
      # ---- |-Obtain the fitting dataframe from data ----
      all_season_fit_wide = wrangle_fit_df(params,data,all_season_country,country_short_input,target_input)
      # ---- |-Make stan list ----
      stan_list = make_stan_list(params,data,all_season_fit_wide,country_short_input,vax_country,pop_country)
      # replace by fake data
      if (F) stan_list = generate_ili_epi_test(par = c(NA),stan_list)
      stan_list$ili_obs_fit %>% sum() # sim:7879760 (nonoise: 8385360), IT:106110
      path_fit = paste0("../Big data/multiseason_age_vax",target_input,country_short_input,".Rdata")
      pr=paste("> Now fitting:",target_input,"for",country_short_input,"... "); cleancat(green(pr))
      
      # ---- |-Fit ----
      if (F) fitout=fit_with_eabc(params,stan_list)
      mod_path='./stan/SIR_multiseason_age_vax.stan'
      if (T) fitout=fit_with_stan(params,stan_list,mod_path=mod_path,all_season_fit_wide)
      
      save(fitout,stan_list,file = path_fit)
      pr=paste("> Fitting:",target_input,"for",country_short_input,"Done \n"); cleancat(green(pr))
      
      modl[[country_short_input]] = fitout
    }
    end_time <- Sys.time() # 
    
    # add stuff to the output list
    df_out$submission = NULL
    df_out$other = modl
  }
  
  if (F) source("code/01_main_supporting/calling_other_models.R") # once you want to call additional models
  
  pr=paste("> Method completed:",round(end_time - start_time,2),"sec \n"); cat(green(pr))
  print(end_time - start_time)
  #### output 
  return(df_out)
}

plot_add_season = function(all_season) {
  p = all_season %>% 
    filter(season%in%c(params$SIR_multiseason$seasons_include) ) %>% 
    unnest(respicompass_ili_plus) %>% 
    ggplot(aes(date,value,color=season)) + geom_line() + 
    facet_wrap(~country_short,scales="free_y") 
  return(p)
}
