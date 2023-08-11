run_flu_models = function(params=NULL, data=NULL){
  
  # Columns in the resulting df_out:
  # model, country_short, date, agegroup, target, value, scenario_tag
  # prediction_type (sample or quantile), sample_or_quantile
  
  df_out = list()
  
  if ( "SIR_simple" %in% params$models_to_run ){ # Old DK model :)
    # Run DK model
    country_short = "AT"
    agegroup = "age_total"
    target = "ILIcases"
    scenario_tag = "A"
    
    df = model_SIR_simple( params, data, country_short, agegroup, target, scenario_tag)
    
    df_out %<>% bind_rows(df) # Add DK model to the df_out
  }
  
  return(df_out)
}