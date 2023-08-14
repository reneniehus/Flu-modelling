last_year_burden = function( params, data, country_short_input, scenario_tag ){
  df_out = NULL
  # create the dataframe for fitting
  data_mock = data %>% 
    filter(country_short == country_short_input, 
           target == params$SIR_simple$target, 
           agegroup == params$SIR_simple$agegroup) # 
  
  data_mock_fit %<>% filter(date>"2022-06-01") # FIXME: remove hardcoding
  # Pcreate the dataframe for projecting
  data_mock_project = data_mock
  data_mock_project$date = data_mock_project$date+365
  
  # loop day by day through last season
  for (d_i in 1:nrow(data_mock_fit)) { 
    # FIXME: more logical to look through project df and then refer back to "mean" year and "sd" years
    date_i = data_mock_fit$date[d_i]
    # use the log-mean of last season
    focus_season_value = data_mock_fit$value[data_mock_fit$date==date_i]
    log_mean = log( focus_season_value+0.0001 ) # FIXME:hardcoding
    
    # get other seasons for variance
    values_vec = rep(0,10) # FIXME: hardcoded
    for (year_i in 1:10){
      date_past = date_i - year_i*365
      date_past_closest = data_mock$date[ which.min(abs(date_past-data_mock$date)) ]
      # FIXME: add a warning if dates are too far
      values_vec[year_i] = data_mock$value[ data_mock$date == date_past_closest ]
    }
    log_sd = sd(log(values_vec))
    
    # sampling
    for (sampl_i in 1:100) {
      sample_value = exp( rnorm(n=1, mean=log_mean , sd=log_sd ) )
      # FIXME: seed setting
      df = tibble(
        model = "last_year_burden",
        country_short = country_short_input,
        agegroup = params$SIR_simple$agegroup,
        target = params$SIR_simple$target,
        scenario_tag = scenario_tag,
        prediction_type = "sample",
        sample_or_quantile = sampl_i,
        value = sample_value
      )
      df_out %<>% bind_rows(df)
    }
    
  }
  #
} 