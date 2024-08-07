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