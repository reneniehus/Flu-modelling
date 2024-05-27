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
    df = model_SIR_simple( params, data$erviss_ili_ari, country_short_input, date_v_fit, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  if ( "SIR_simple_through_data" %in% params$models_to_run ){
    df_collect = list()
    df_i = 1
    country_short_input_v = unique(data$erviss_ili_ari$country_short) 
    
    for (country_short_input_i in country_short_input_v) {
      start_year = data$erviss_ili_ari %>% filter(country_short==country_short_input_i) %>% pull(date) %>% min() %>% year() %>% as.numeric()
      while(start_year<=2023) {
        season = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,"-07-01"))
        end_date = ymd(paste0(start_year+1,"-05-01"))
        start_year = start_year +1 
        date_v_fit = seq(from=start_date,to=end_date,by="day")
        
        # test filtering
        # ili
        data$erviss_ili_ari %>% 
          filter(country_short == country_short_input_i, 
                 target == params$SIR_simple$target, 
                 agegroup == params$SIR_simple$agegroup) %>% 
          filter( date%in%date_v_fit ) -> xfit
        # detections
        data$erviss_detect_sent %>% 
          filter(country_short == country_short_input_i,date%in%date_v_fit) %>% 
          filter(pathogen=="Influenza") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age) -> xdet_sent
        ntests = xdet_sent %>% filter(indicator=="tests") %>% summarise(msum=sum(value)) %>% pull(msum)
        # detections
        data$erviss_detect_nonsent %>% 
          filter(country_short == country_short_input_i,date%in%date_v_fit) %>%
          filter(pathogen=="Influenza") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age) -> xdet_nonsent
        ntests_nonsent = xdet_nonsent%>% filter(indicator=="tests")%>% summarise(msum=sum(value)) %>% pull(msum)
        
        xfit %>% ggplot(aes(date,value))+geom_line()
        xdet_sent %>% ggplot(aes(date,value))+geom_line()
        # if ( nrow(xfit) < 10 ) next;
        # sum_inc = sum(xfit$value) ; if ( sum_inc < 300 ) next;
        pr=paste("> Running:",country_short_input_i,"| season:",season,
                 "| sum inc:",sum_inc,"\n"); cat(green(pr))
        
        # put data together
        df_collect[[df_i]] = tibble(
          country_short=country_short_input_i,
          season=season,
          ili_sum= sum(xfit$value),
          sentinel_tests=ntests,
          nonsentinel_test=ntests_nonsent
          )
        #
        df_i = df_i + 1
      }
    }
    all_season = bind_rows(df_collect)
    all_season %>% mutate(country=EU_long(country_short)) %>%  
      group_by(country) %>% 
      summarise(ILI=mean(ili_sum) %>% round(),
                ILI_min=min(ili_sum) %>% round(),
                sent_tests=mean(sentinel_tests) %>% round(),
                sent_tests_min=min(sent_tests) %>% round(),
                nonsent_test=mean(nonsentinel_test) %>% round(),
                nonsent_test_min=min(nonsentinel_test) %>% round()) %>% View()
  }
  
  
  if ( "SIR_simple_r0_variation" %in% params$models_to_run ){ # Old DK model :)
    
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
          filter( date%in%date_v_fit ) -> xfit
        xfit %>% ggplot(aes(date,value))+geom_line()
        if ( nrow(xfit) < 10 ) next;
        sum_inc = sum(xfit$value) ; if ( sum_inc < 300 ) next;
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
