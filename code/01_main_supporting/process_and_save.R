process_and_save = function(params=NULL, models_out=NULL){
  
  # file naming (https://docs.google.com/document/d/13adcxpPdlDvJM5eiFSkMzlWMTcwsx6lVjY25JA26iS4/edit)
  # season_cycle_round_id>-<team>-<model>.parquet (Ex. 2024_2025_1_FLU2-ISI-GLEAM.parquet)

  write.fst(models_out, path = params$path_save_results)
  
  # Basic plot(s):
  # Show number of samples per model and date 
  avail_forecasts = 
    avail_forecasts(models_out %>% mutate(prediction=value, true_value=sample_or_quantile, sample=sample_or_quantile),
                    by=c("scenario_tag","country_short","date","model"),
                    collapse = "sample_or_quantile")
  p = plot_avail_forecasts(
    avail_forecasts, x="date", show_numbers = F
  ) + facet_grid("scenario_tag"~"country_short")
  ggsave(filename = paste0(params$path_save_figures,"sample_overview.jpg"), plot = p)
}