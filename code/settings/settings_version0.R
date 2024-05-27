settings = function() {
  params = list()
  
  # ---- |-Disease parameters ----
  params$rate_infectious = 0.2777778
  params$Rnull = 2.0; # https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  # ---- |-Folder paths ----
  params$path_save_results = "./output/output.fst"
  params$path_save_figures = "./output/figures/"
  
  # ---- |-Data source ----
  
  # ---- |-Countries ----
  
  # ---- |-Model-specific  settings ----
  params$models_to_run = c("SIR_simple_r0_variation") # "SIR_simple"
  
  # Settings for SIR_simple
  params$SIR_simple$target = "ILIconsultationrate"
  params$SIR_simple$agegroup = "age_total"
  # Settings for last_year_burden
  params$last_year_burden$target = "ILIconsultationrate"
  params$last_year_burden$agegroup = "age_total"
  # Settings for arima simple
  params$arima_simple$target = "ILIconsultationrate"
  params$arima_simple$agegroup = "age_total"
  
  # ---- |-Fitting and uncertainty ----
  
  # ---- |-Flu scenarios ----
  params$scenarios = c("baseline")
  params$scenarios = c("baseline")
  
  return(params)
}

 



