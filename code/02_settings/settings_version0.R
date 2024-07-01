settings = function() {
  params = list()
  
  params$debug = T
  # ---- |-Names/identifiers ----
  params$scenario_round_id = "2024_2025_1_FLU1"
  params$scenario_team = "ECDC"
  # ---- |-Disease parameters ----
  params$rate_infectious = 0.2777778
  params$Rnull = 2.0 # https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  # vaccine parameters
  params$ve_spread = 0.20 # assuming a similar size of effect as for ve_inf
  params$ve_inf = 0.20 #
  params$ve_ili_cond_inf = 0.25 # 
  # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf) 
  params$ve_ili = 1-(1-params$ve_inf)*(1-params$ve_ili_cond_inf) # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf) 
  
  params$ve_inf = 0.3 # vaccine efficacy on infectiousness
  params$ve_susc = 0.3 # vaccine efficacy on susceptability
  params$ve_severe = 0.6 # vaccine efficacy on severity, given infection
  
  # ---- |-Data ----
  params$latest_start_year = 2023 # if the last full season is 2023/24, put 2023
  params$season_start_monthday = "-07-01"
  params$season_end_monthday = "-06-30"
  
  # ---- |-Countries ----
  
  # ---- |-Model-specific  settings ----
  params$models_to_run = c("SIR_simple_multi_season") # "SIR_simple","SIR_simple_r0_variation","SIR_simple_multi_season
  
  # Settings for SIR_simple
  params$SIR_simple$target = "ILIconsultationrate"
  params$SIR_simple$agegroup = "age_total"
  # Settings for SIR_multiseason
  params$SIR_simple_multi_season$target = c("respicompass_ili_plus") # c("ili","ili_typing_sentinel","ili_typing_all","respicompass_ili_plus")
  params$SIR_multiseason$seasons_exclude = c("2019/2020","2020/2021","2021/2022") # those impacted by COVID-19 acute phase
  params$SIR_multiseason$seasons_include = c("2017/2018","2018/2019","2023/2024") # 2017-2018, 2018-2019, and 2023-2024
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
  
  # ---- |-Folder paths ----
  # season_cycle_round_id>-<team>-<model>.parquet (Ex. 2024_2025_1_FLU2-ISI-GLEAM.parquet)
  for (model_id in params$models_to_run) params$path_save_results[[model_id]] = paste0("./output/",
                                                                                       params$scenario_round_id,"-",
                                                                                       params$scenario_team,"-",
                                                                                       model_id,".parquet")
  
  params$path_save_figures = "./output/figures/"
  
  return(params)
}

 



