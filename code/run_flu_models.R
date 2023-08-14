run_flu_models = function(params=NULL, data=NULL){
  
  # Columns in the resulting df_out:
  # model, country_short, date, agegroup, target, value, scenario_tag
  # prediction_type (sample or quantile), sample_or_quantile
  
  df_out = list()
  
  if ( "SIR_simple" %in% params$models_to_run ){ # Old DK model :)
    # Run SIR_simple model
    country_short = "AT"
    scenario_tag = "A"
    
    df = model_SIR_simple( params, data, country_short, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  
  
  if ( "last_year_burden" %in% params$models_to_run ){ # O
    # un last_year_burden model
    country_short = "AT"
    scenario_tag = "A"
    
    df = last_year_burden( params, data, country_short, scenario_tag)
    last_year_burden = function( params, data, country_short, scenario_tag ){
      # create the dataframe for fitting
      data_mock = data %>% 
        filter(country_short == country_short_input, 
               target == params$SIR_simple$target, 
               agegroup == params$SIR_simple$agegroup) # 
      data_mock %<>% filter(date>"2022-06-01")
      data_mock_fit = data_mock
      # Pcreate the dataframe for projecting
      data_mock_project = data_mock
      data_mock_project$date = data_mock_project$date+365
      
      # loop day by day through last season
      for (d_i in 1:nrow(data_mock_fit)) {
        # use the log-mean of last season
        # use the log_se of other seasons
      }
      #
    }
    
    
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  return(df_out)
}

#########################

# Create an output dataframe
df_out = fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>%
  ungroup() %>%
  left_join(data_mock %>% select(date) %>% mutate(t_vw = 1:n()), by="t_vw") %>%
  mutate(model = "SIR_simple",
         country_short = country_short,
         agegroup = agegroup,
         target = target,
         scenario_tag = scenario_tag,
         prediction_type = "sample"
  ) %>%
  rename(value = .value,
         sample_or_quantile = .draw) %>%
  select(-.chain, -.iteration, -.variable, -t_vw) 