run_flu_models = function( params=NULL , data=NULL ){
  
  # ---- |-initiating desired output list ----
  df_out = list(
    time_of_execution = now(),    # time-stamp
    df_for_submission = NULL,     # for each model a clean dataframe following submission format 
    output_other = NULL           # for each model additional output 
  )
  
  # ---- |-Run analysis ----
  if ( "SIR_simple_r0_variation" %in% params$models_to_run ){ 
    
    # data
    all_season = data_into_all_season(data,params,withforce=F)
    dat = data$epi$erviss_ili_ari
    country_short_input_v = unique(dat$country_short) # country_short_input_v = "AT" # for quick run
    scenario_tag = "A"
    target_input = "ili_typing_sentinel"
    
    # prepare model input
    df_collect = list()
    df_i = 1
    
    start_time <- Sys.time()
    
    for (country_short_input_i in country_short_input_v) { # country_short_input_i = country_short_input_v[1]
      
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
    } # country loop
    end_time <- Sys.time()
    (end_time - start_time) # 1.6 hours
    
    df_collect %>% bind_rows -> x
    if (params$debug==F) write_csv(x,file="output/rt_season_country.csv")
  }
  
  return(df_out)
}
