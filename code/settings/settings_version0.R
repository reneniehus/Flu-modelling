settings = function() {
  params = list()
  # ---- |-Flu parameters ----
  params$rate_infectious = 0.2777778
  params$Rnull = 2.0; # https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  
  params$path_save_results = "./output/output.fst"
  params$path_save_figures = "./output/figures/"
  # ---- |-Data source ----
  
  # ---- |-Countries ----
  
  # ---- |-Models ----
  params$models_to_run = c("SIR_simple") # "SIR_simple",
  # Settings for SIR_simple
  params$SIR_simple$target = "ILIcases"
  params$SIR_simple$agegroup = "age_total"
  # Settings for last_year_burden
  params$last_year_burden$target = "ILIcases"
  params$last_year_burden$agegroup = "age_total"
  
  params$arima_simple$target = "ILIcases"
  params$arima_simple$agegroup = "age_total"
  
  # ---- |-Fitting and uncertainty ----
  
  # ---- |-Simulation setting ---- 
  
  # ---- |-Flu scenerios ----
  return(params)
}

 



