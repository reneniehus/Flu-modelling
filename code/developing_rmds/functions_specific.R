get_weekly_daygroups = function(df_Rt_value_pop){
  df_out = list()
  n = nrow(df_Rt_value_pop)
  
  n_full_weeks = n - n %% 7
  nweeks = n_full_weeks/7
  df_out$nweeks = (nweeks-1)
  d_id = c(1:n_full_weeks)
  week_day_map = rep(c(1:nweeks),each=7)
  
  weekly_daygroups_t = list()
  weekly_daygroups_tlast = list()
  for (i in 1:(nweeks-1) ) {
    weekly_daygroups_t[[i]] = d_id[week_day_map==(i+1) ]
    weekly_daygroups_tlast[[i]] = d_id[week_day_map==(i) ]
  }
  df_out$n_full_weeks = n_full_weeks
  df_out$weekly_daygroups_t = weekly_daygroups_t
  df_out$weekly_daygroups_tlast = weekly_daygroups_tlast
  
  return(df_out)
}

weekly_make_daily = function(df_ii_weekly){
  size_of_week = 7
  adding_half_week_at_tail = 3
  date_range = range(df_ii_weekly$date)
  daily_v = seq(from=date_range[1],to=date_range[2]+adding_half_week_at_tail,by="day")
  df_out = tibble(date=daily_v) %>% left_join(df_ii_weekly,by = join_by(date))
  df_out = df_out %>% 
    mutate(value_approx = na.approx(value,na.rm = F)/size_of_week ) %>% 
    fill(value,.direction = "down") %>% fill(location,.direction = "down") %>% 
    mutate( value=rollmean(value,k=size_of_week,fill=NA,align="right") ) %>% 
    mutate( value=value/size_of_week) 
  if (F) {
    df_out %>%  ggplot(aes(x=date,y=value)) + geom_line() + 
      geom_line(aes(y=value_approx),col="red") + scale_y_log10()
  } 
  return(df_out)
}