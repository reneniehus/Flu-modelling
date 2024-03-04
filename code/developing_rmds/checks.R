# check df_ii
no_duplication_no_gaps_right_order = function(df_ii){
  #
  any_duplicated = (df_ii$date) %>% duplicated() %>% any()
  if (any_duplicated==T) {pr=paste("df_ii has duplications in date","!!!\n"); cat(red(pr))}
  #
  num_unique_gaps = (df_ii$date) %>% diff() %>% unique() %>% length()
  if (num_unique_gaps!=1) {pr=paste("df_ii has uneven gaps","!!!\n"); cat(red(pr))}
  #
  sorted_already =  all.equal(df_ii$date,sort(df_ii$date))
  if (sorted_already==F) {pr=paste("df_ii not sorted correctly","!!!\n"); cat(red(pr))}
  #
  any_na_value = (df_ii$value) %>% is.na() %>% any()
  if (any_na_value==T) {pr=paste("df_ii contains NAs in value column","!!!\n"); cat(red(pr))}
  #
  if (any_duplicated==F & num_unique_gaps==1 & sorted_already==T & any_na_value==F) {
    pr=paste("df_ii has NO duplications, has evenly spaced gaps, is arranged correctly, without NAs in value",":-) \n"); cat(green(pr))
  }
}

# check df_transmission

check_df_transmission = function(df_transmission){
  Rt_null = is.null(df_transmission$Rt)
  if (Rt_null==T) {pr=paste("df_transmission has no Rt column yet","!!!\n"); cat(red(pr))}
  
  date_diffs = df_transmission$date %>% diff() %>% unique()
  if (date_diffs!=1) {pr=paste("df_transmission is not daily","!!!\n"); cat(red(pr))}
  
  n_na = is.na(df_transmission$Rt) %>% sum()
  n_all = length(df_transmission$Rt)
  if (Rt_null==F) {pr=paste("df_transmission contains Rt with",n_na,"(",round(n_na/n_all*100),"%)","NA values!!!\n"); cat(green(pr))}
}

# check df_Rt_value_pop
not_too_short_or_long = function(df_Rt_value_pop){
  date_range = df_Rt_value_pop$date %>% range()
  date_range_in_days = as.numeric(diff(date_range))
  very_long = (date_range_in_days - 365) > 2*30
  if (very_long==T) {pr=paste("not_too_short_or_long is very long, check","!!!\n"); cat(red(pr))}
  very_short = (date_range_in_days - 365) < -2*30
  if (very_short==T) {pr=paste("not_too_short_or_long is very short, check","!!!\n"); cat(red(pr))}
  
  if (very_long==F & very_short==F ) {
    pr=paste("df_Rt_value_pop has a good length",":-) \n"); cat(green(pr))
  }
  
}

# check stan_list
checking_stan_list = function(stan_list){
  #
  pr=paste("Model will be fit to",stan_list$n_location,"location:\n"); cat(blue(pr))
  pr=paste(stan_list$location_id_raw$value); cat(white(pr))
  #
  pr=paste("Model will be fit using",stan_list$n,"data points.\n"); cat(blue(pr))
  # priors
  if (F) {pr=paste("You prior for X should be closer to",stan_list$n,".\n"); cat(red(pr)) }
}