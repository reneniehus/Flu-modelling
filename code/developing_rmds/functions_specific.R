get_weekly_daygroups = function(df_all){
  df_out = list()
  n = nrow(df_all)
  
  n_full_weeks = n - n %% 7
  nweeks = n_full_weeks/7
  df_out$nweeks = (nweeks-1)
  d_id = c(1:n_full_weeks)
  week_day_map = rep(c(1:weeks),each=7)
  
  weekly_daygroups_t = list()
  weekly_daygroups_tlast = list()
  for (i in 1:(nweeks-1) ) {
    weekly_daygroups_t[[i]] = d_id[week_day_map==(i+1) ]
    weekly_daygroups_tlast[[i]] = d_id[week_day_map==(i) ]
  }
  df_out$weekly_daygroups_t = weekly_daygroups_t
  df_out$weekly_daygroups_tlast = weekly_daygroups_tlast
  
  return(df_out)
}
