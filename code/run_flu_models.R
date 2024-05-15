run_flu_models = function(params=NULL, data=NULL){
  
  # Columns in the resulting df_out:
  # model, country_short, date, agegroup, target, value, scenario_tag
  # prediction_type (sample or quantile), sample_or_quantile
  
  df_out = list()
  
  if ( "SIR_simple" %in% params$models_to_run ){ # Old DK model :)
    # prepare for model
    country_short_input = "AT"
    scenario_tag = "A"
    date_v_fit = seq(from=ymd("2022-10-05"),to=ymd("2023-05-01"),by="day")
    # Run SIR_simple model
    df = model_SIR_simple( params, data, country_short_input, date_v_fit, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  
  if ( "SIR_simple_r0_variation" %in% params$models_to_run ){ # Old DK model :)
    
    # prepare for model
    df_collect = list()
    df_i = 1
    scenario_tag = "A"
    country_short_input_v = unique(data$country_short) # country_short_input_v = "AT" # for quick run
    for (country_short_input_i in country_short_input_v) {
      country_short_input = country_short_input_i
      
      start_year = data %>% filter(country_short==country_short_input) %>% pull(date) %>% min() %>% year() %>% as.numeric()
      while(start_year<=2022) {
        season = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,"-08-01"))
        end_date = ymd(paste0(start_year+1,"-05-01"))
        start_year = start_year +1 
        date_v_fit = seq(from=start_date,to=end_date,by="day")
        
        # test filtering
        data %>% 
          filter(country_short == country_short_input, 
                 target == params$SIR_simple$target, 
                 agegroup == params$SIR_simple$agegroup) %>% 
          filter( date%in%date_v_fit ) -> xfit
        xfit %>% ggplot(aes(date,value))+geom_line()
        if ( nrow(xfit) < 10 ) next;
        sum_inc = sum(xfit$value) ; if ( sum_inc < 300 ) next;
        pr=paste("> Running:",country_short_input,"| season:",season,"| sum inc:",sum_inc,"\n"); cat(green(pr))
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
  
  if ( "last_year_burden" %in% params$models_to_run ){ # O
    # prepare for run
    country_short_input = "AT"
    scenario_tag = "A"
    # run last_year_burden model
    df = last_year_burden( params, data, country_short_input, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  if ( "arima_simple" %in% params$models_to_run ){ # O
    # prepare for run
    country_short_input = "AT"
    scenario_tag = "A"
    # run last_year_burden model
    df = arima_simple( params, data, country_short_input, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  return(df_out)
}
