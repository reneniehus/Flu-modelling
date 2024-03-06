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
