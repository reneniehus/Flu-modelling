run_flu_models = function(params=NULL, data=NULL){
  
  # Columns in the resulting df_out:
  # model, country_short, date, agegroup, target, value, scenario_tag
  # prediction_type (sample or quantile), sample_or_quantile
  
  df_out = list()
  
  if ( "SIR_simple" %in% params$models_to_run ){ # Old DK model :)
    # prepare for model
    country_short_input = "AT"
    scenario_tag = "A"
    # Run SIR_simple model
    df = model_SIR_simple( params, data, country_short_input, scenario_tag)
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  
  
  if ( "last_year_burden" %in% params$models_to_run ){ # O
    # prepare for run
    country_short_input = "AT"
    scenario_tag = "A"
    
    # run last_year_burden model
    df = last_year_burden( params, data, country_short_input, scenario_tag)
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