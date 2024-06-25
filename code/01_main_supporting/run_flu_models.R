run_flu_models = function( params=NULL , data=NULL ){
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Notes on output requirements ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # Required columns in df_for_submission (as per: https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit):
  # model_id
  # round_id ["2024_2025_1_FLU1"]
  # scenario_id ["A","B"], target [allowed targets], location ["DE","FR"] 
  # pop_group ["0-12","13-65"], horizon [week integer], target_end_date [Date string ('YYYY-MM-DD')]
  # output_type ["sample"], output_type_id [string: "1","2","3",...], value [float limited to 2 decimals]
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Initiating desired output list ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_out = list(
    time_of_execution = now(),    # time-stamp
    df_for_submission = NULL,     # for each model a clean dataframe following submission format (see above)
    output_other = NULL           # for each model additional output 
  )
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Running selected models ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  if ( "SIR_simple" %in% params$models_to_run ){ 
    # prepare for model
    country_short_input = "AT"
    date_v_fit = seq(from=ymd("2022-10-05"),to=ymd("2023-05-01"),by="day")
    # Run SIR_simple model
    df = model_SIR_simple( params, dat=data$epi$erviss_ili_ari, country_short_input, date_v_fit )
    # add data
    df_out$df_for_submission[["SIR_simple"]] = df
    df_out$output_other[["SIR_simple"]] = list(
      date_v_fit = date_v_fit
    )
  }
  
  if ( "SIR_simple_multi_season" %in% params$models_to_run ){
    
    pr=paste("Initiating SIR_simple_multi_season \n"); cat(green(pr))
    
    # ---- |-Data ----
    all_season = data_into_all_season(data,params,withforce=F)
    country_short_input_v = all_season %>% filter_log(ili_sum>0) %>% pull(country_short) %>% unique()
    target_input_v = params$SIR_simple_multi_season$target
    scenario_tag = "A"
    
    modl = list()
    start_time <- Sys.time()
    
    # ---- |-Run model ----
    for ( target_input in target_input_v ) { # target_input=target_input_v[1]
      for (country_short_input in country_short_input_v ) { # country_short_input=country_short_input_v[1]
        # population
        pop_country = data$demography$population_pyramid %>% 
          filter(country==country_short_input) %>% pull(population) %>% sum()
        if (country_short_input=="GR") pop_country = 10.43*1e6
        vax_country = data$vax$data_vax %>% filter( location_name == EU_long(country_short_input) ) # vaccination data for a country
        if (nrow(vax_country) != 1) stop("Vaccination data is wrong format: either no data or too many rows")
        # run model
        modl[[target_input]][[country_short_input]] = model_SIR_multiseason( params , 
                                                                             all_season=all_season , 
                                                                             target_input, 
                                                                             country_short_input,
                                                                             pop_country,
                                                                             vax_country)
      }
    }
    end_time <- Sys.time() # 5 hrs
    pr=paste("> Method run:",round(end_time - start_time,2),"sec \n"); cat(green(pr))
    save(modl,file = "../Big data/modl.Rdata")
    
    # ---- |-Analyse model output ----
    mdf_all = NULL
    for (target_input in target_input_v)
    for (country_short_input in country_short_input_v){
      fit = modl[[target_input]][[country_short_input]]$fit
      stan_list = modl[[target_input]][[country_short_input]]$stan_list
      season_df = stan_list$season_id_raw %>% 
        rename(season_numeric=name,season=value)
      mdf = crossing( target=target_input,country=country_short_input , 
                      season_df )
      mdf$prop_severe = precis(fit,"prop_severe",depth=2)$mean
      
      sir_names = precis(fit,"SIR_ini",depth=3) %>% rownames()
      pat <- "(\\d)+"
      season_id = as.numeric(str_extract(sir_names, pat))
      sir_id = rep(c("S","I","R"),nrow(mdf))
      mdf$S_ini = precis(fit,"SIR_ini",depth=3)$mean[sir_id=="S"]
      mdf$R_ini = precis(fit,"SIR_ini",depth=3)$mean[sir_id=="R"]
      mdf$reciprocal_phi = precis(fit,"reciprocal_phi",depth=1)[1,1]
      
      mdf$prop_severe_rhat = precis(fit,"prop_severe",depth=2)$Rhat4
      mdf$S_ini_rhat = precis(fit,"SIR_ini",depth=3)$Rhat4[sir_id=="S"]
      mdf$R_ini_rhat = precis(fit,"SIR_ini",depth=3)$Rhat4[sir_id=="R"]
      mdf$reciprocal_phi_rhat = precis(fit,"reciprocal_phi",depth=1)[6,1]
      mdf_all = bind_rows(mdf_all,mdf)
    }
  
    mdf_all %>% filter(!target=="ili_typing_all") %>% 
      ggplot(aes(season,logit(prop_severe),fill=target)) + geom_boxplot()
    
    mdf_all %>% filter(!target=="ili_typing_all") %>% 
      ggplot(aes(season,logit(prop_severe),group=country)) + geom_line() + 
      facet_wrap(~target)
    
    mdf_all %>% filter(!target=="ili_typing_all") %>% 
      ggplot(aes(season,logit(S_ini),fill=target)) + geom_boxplot()
    mdf_all %>% filter(!target=="ili_typing_all") %>% 
      ggplot(aes(season,logit(R_ini),fill=target)) + geom_boxplot()
    
    mdf_all %>% 
      ggplot(aes(target,(reciprocal_phi))) + geom_boxplot()
    
    mcountry ="AT"
    p1=modl[["ili"]][[mcountry]]$pdata
    p2=modl[["ili_typing_sentinel"]][[mcountry]]$pdata
    p3=modl[["ili_typing_all"]][[mcountry]]$pdata
    p1/p2/p3
    
    df = NULL
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  if ( "SIR_simple_r0_variation" %in% params$models_to_run ){ 
    
    # data
    all_season = data_into_all_season(data,params,withforce=F)
    
    # prepare model input
    df_collect = list()
    df_i = 1
    scenario_tag = "A"
    target_input = "ili_typing_sentinel"
    
    country_short_input_v = unique(dat$country_short) # country_short_input_v = "AT" # for quick run
    start_time <- Sys.time()
    for (country_short_input_i in country_short_input_v) {
      
      pop_country = data$demography$population_pyramid %>% 
        filter(country==country_short_input_i) %>% pull(population) %>% sum()
      if (country_short_input_i=="GR") pop_country = 10.43*1e6
      
      start_year = dat %>% filter(country_short==country_short_input_i) %>% pull(date) %>% min() %>% year() %>% as.numeric()
      while(start_year<=2022) {
        season = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,"-07-01"))
        end_date = ymd(paste0(start_year+1,"-05-01"))
        start_year = start_year +1 
        date_v_fit = seq(from=start_date,to=end_date,by="day")
        
        # test filtering
        dat %>% 
          filter(country_short == country_short_input_i, 
                 target == params$SIR_simple$target, 
                 agegroup == params$SIR_simple$agegroup) %>% 
          filter( date%in%date_v_fit ) -> xinc_iliari
        xinc_iliari %>% ggplot(aes(date,value))+geom_line()
        if ( nrow(xinc_iliari) < 10 ) next;
        sum_inc = sum(xinc_iliari$value) ; if ( sum_inc < 300 ) next;
        pr=paste("> Running:",country_short_input_i,"| season:",season,"| sum inc:",sum_inc,"\n"); cat(green(pr))
        df_collect[[df_i]] = model_SIR_simple_r0( params, all_season=all_season , target_input, pop_country, country_short_input=country_short_input_i, date_v_fit,season )
        df_i = df_i + 1
        
      } # season loop
    } # 
    end_time <- Sys.time()
    (end_time - start_time) # 1.6 hours
    
    
    if (T){
      df_collect %>% bind_rows -> x
      write_csv(x,file="code/03_special_analyses/rt_season_country.csv")
    }
    x = read_csv(file="code/03_special_analyses/rt_season_country.csv")
    x = df_collect %>% bind_rows()
    (x$Rnull_eff) %>% min()
    rnull_mu = x$Rnull_eff %>% median()
    rnull_quant = x$Rnull_eff %>% quantile(probs=c(0.2,0.8))
    ((rnull_quant/rnull_mu )-1)*100 # -11.960425   9.747009 
  }
  
  if ( "last_year_burden" %in% params$models_to_run ){ # 
    # prepare for run
    country_short_input = "AT"
    scenario_tag = "A"
    # run last_year_burden model
    df = last_year_burden( params, data, country_short_input, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  if ( "arima_simple" %in% params$models_to_run ){ # 
    # prepare for run
    country_short_input = "AT"
    scenario_tag = "A"
    # run last_year_burden model
    df = arima_simple( params, data, country_short_input, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  return(df_out)
}
