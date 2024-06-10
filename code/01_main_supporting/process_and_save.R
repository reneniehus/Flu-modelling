process_and_save = function(params=NULL, models_out=NULL){
  
  # file naming (https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit)
  # season_cycle_round_id>-<team>-<model>.parquet (Ex. 2024_2025_1_FLU2-ISI-GLEAM.parquet)
  for (mod_id in params$models_to_run){
    write.fst(models_out$df_for_submission[[mod_id]], path = params$path_save_results[[mod_id]])
  }
  
  # Basic plot(s):
  # Show number of samples per model and date 
  mdf = NULL
  for (mod_id in params$models_to_run){
    mdf = bind_rows( mdf,models_out$df_for_submission[[mod_id]] )
  }
   
  avail_forecasts = 
    avail_forecasts(mdf %>% mutate(prediction=value, true_value=output_type_id, sample=output_type_id,model=model_id),
                    by=c("scenario_id","location","target_end_date","model"),
                    collapse = "sample_or_quantile")
  p = plot_avail_forecasts(
    avail_forecasts, x="target_end_date", show_numbers = F
  ) + facet_grid("scenario_tag"~"country_short")
  
  ggsave(filename = paste0(params$path_save_figures,"sample_overview.jpg"), plot = p)
  
}