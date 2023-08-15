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
