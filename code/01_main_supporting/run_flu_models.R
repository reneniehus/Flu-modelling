run_flu_models = function(params=NULL, data=NULL){
  
  # Columns in the resulting df_out (as per: https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit):
  # model_id
  # round_id ["2024_2025_1_FLU1"]
  # scenario_id ["A","B"], target [allowed targets], location ["DE","FR"] 
  # pop_group ["0-12","13-65"], horizon [week integer], target_end_date [Date string ('YYYY-MM-DD')]
  # output_type ["sample"], output_type_id [string: "1","2","3",...], value [float limited to 2 decimals]
  
  df_out = list()
  
  if ( "SIR_simple" %in% params$models_to_run ){ 
    # prepare for model
    country_short_input = "AT"
    scenario_tag = "A"
    date_v_fit = seq(from=ymd("2022-10-05"),to=ymd("2023-05-01"),by="day")
    # Run SIR_simple model
    df = model_SIR_simple( params, data$erviss_ili_ari, country_short_input, date_v_fit, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  if ( "SIR_simple_multi_season" %in% params$models_to_run ){
    # prepare data
    all_season = data_into_all_season(data,params,withforce=F)
    
    modl = list()
    scenario_tag = "A"
    target_input_v = c("ili","ili_typing_sentinel","ili_typing_all")
    country_short_input_v = all_season %>% filter_log(ili_sum>0) %>% pull(country_short) %>% unique()
    for (target_input in target_input_v) {
      for (country_short_input in country_short_input_v ) {
        modl[[target_input]][[country_short_input]] = model_SIR_multiseason( params , all_season=all_season , target_input, country_short_input, scenario_tag)
      }
    }
    save(modl,file = "../Big data/modl.Rdata")
    
    
    mcountry ="AT"
    p1=modl[["ili"]][[mcountry]]$pdata
    p2=modl[["ili_typing_sentinel"]][[mcountry]]$pdata
    p3=modl[["ili_typing_all"]][[mcountry]]$pdata
    p1/p2/p3
    
    df = NULL
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  
  if ( "SIR_simple_r0_variation" %in% params$models_to_run ){ 
    # prepare for model
    df_collect = list()
    df_i = 1
    scenario_tag = "A"
    country_short_input_v = unique(data$erviss_ili_ari$country_short) # country_short_input_v = "AT" # for quick run
    for (country_short_input_i in country_short_input_v) {
      country_short_input = country_short_input_i
      
      start_year = data %>% filter(country_short==country_short_input_i) %>% pull(date) %>% min() %>% year() %>% as.numeric()
      while(start_year<=2022) {
        season = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,"-07-01"))
        end_date = ymd(paste0(start_year+1,"-05-01"))
        start_year = start_year +1 
        date_v_fit = seq(from=start_date,to=end_date,by="day")
        
        # test filtering
        data %>% 
          filter(country_short == country_short_input_i, 
                 target == params$SIR_simple$target, 
                 agegroup == params$SIR_simple$agegroup) %>% 
          filter( date%in%date_v_fit ) -> xinc_iliari
        xinc_iliari %>% ggplot(aes(date,value))+geom_line()
        if ( nrow(xinc_iliari) < 10 ) next;
        sum_inc = sum(xinc_iliari$value) ; if ( sum_inc < 300 ) next;
        pr=paste("> Running:",country_short_input_i,"| season:",season,"| sum inc:",sum_inc,"\n"); cat(green(pr))
        df_collect[[df_i]] = model_SIR_simple_r0( params, data, country_short_input, date_v_fit,season, scenario_tag)
        df_i = df_i + 1
      }
    }
    
    if (F){
      df_collect %>% bind_rows -> x
      write_csv(x,file="../Big data/Rt_country_season.csv")
    }
    x = read_csv(file="../Big data/Rt_country_season.csv")
    x = df_collect %>% bind_rows()
    (x$Rnull) %>% min()
    rnull_mu = x$Rnull %>% median()
    rnull_quant = x$Rnull %>% quantile(probs=c(0.2,0.8))
    ((rnull_quant/rnull_mu )-1)*100
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
  
  mout = list(
    df_out,
    multiseason=modl
  )
  
  return(mout)
}
