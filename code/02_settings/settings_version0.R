settings = function() {
  params = list()
  
  # ---- |-Run modes ----
  params$save_submission = F # T: saves the file ready for respicompass, F; will be faster
  
  # debug/fast modes
  params$rapid_stan_fit = F # T: runs scripts with settings that reduce run-time
  params$load_earlyfit  = F # T: saved fits will be used
  
  # ---- |-Resport setting ----
  params$send_report = T
  params$report_recipients = c('rene.niehus@ecdc.europa.eu', 'rene7niehus@gmail.com','rok.grah@ecdc.europa.eu')
  #params$report_recipients = c('rene.niehus@ecdc.europa.eu')
  
  # ---- |-Names/identifiers ----
  params$scenario_round_id = "2024_2025_1_FLU"
  params$scenario_team = "ECDC"
  params$scenario_model = "flumod"
  
  params$four_age_groups = c("0-4","5-14","15-64","65+") # the order is important
  params$proj_start_year = 2024 # when projecting the season 2024/25, put 2024
  
  # ---- |-Disease parameters ----
  params$rate_infectious = 0.2777778
  params$Rnull = 2.0 # https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  
  # immunity parameters
  params$ve_spread = 0.20 # vaccine effect on onward spread when vaccinated individual is infected
  params$ve_inf = 0.25 # vaccine efficacy on becoming infected given exposure
  params$ve_ili_cond_inf = 0.20 # vaccine effect on ILI development given infections
  # vaccine effect on ILI given exposure is the combined effect of ve_inf and ve_ili_cond_inf
  # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf) 
  params$ve_ili = 1-(1-params$ve_inf)*(1-params$ve_ili_cond_inf) # (1-ve_ili) = (1-ve_inf)*(1-ve_ili_cond_inf) 
  
  # ---- |-Data ----
  params$latest_start_year = 2024 # if the last partly/fully observed season is 2024/25, put 2024
  params$season_start_monthday = "-08-01" # initial date of for SIR initiation
  params$season_end_monthday = "-07-31" # end date of SIR process
  
  params$ili_plus_sentinel = c("AT","BE","CZ","DK","EE","FI","FR","IE","IT","NL","NO","PL","SI") # for which countries is ili-plus computation based on sentinel pathogen testing
  params$ili_plus_nonsentinel = c("HR","IS","LV","MT","RO") # for which countries is ili-plus computation based on non-sentinel pathogen testing
  # summer low-activity (where we assume that ILI activity = 0, where NA was reported)
  low_start = "-06-01"
  low_stop  = "-09-01"
  date_v = NULL
  for (year_i in 2010:2030){
    date_v = c(date_v,seq( from=paste0(year_i,low_start) %>% ymd(), to=paste0(year_i,low_stop) %>% ymd(), by="day" ))
  }
  params$summer_low_dates = date_v %>% as_date() 
  
  # ---- |-Simulations ----
  params$simulation_seed = 12
  
  # ---- |-Countries ----
  params$run_countries = "IT"
  
  # ---- |-Model-specific  settings ----
  params$models_to_run = c("SIR_simple_multi_season") # "SIR_simple","SIR_simple_r0_variation","SIR_simple_multi_season
  # Settings for SIR_simple
  params$SIR_simple$target = "ILIconsultationrate"
  params$SIR_simple$agegroup = "age_total"
  # Settings for SIR_multiseason
  params$SIR_simple_multi_season$target = c("erviss_ili_plus") # c("ili","ili_typing_sentinel","ili_typing_all","respicompass_ili_plus","erviss_ili_plus")
  params$SIR_multiseason$seasons_exclude = c("2019/2020","2020/2021","2021/2022") # those impacted by COVID-19 acute phase
  params$SIR_multiseason$seasons_include = c("2017/2018","2018/2019","2023/2024","2024/2025") # 2017-2018, 2018-2019, 2023-2024, and 2024-2025
  params$SIR_multiseason$seasons_baseline = c("2017/2018","2018/2019","2023/2024")
  params$SIR_multiseason$age_groups = c("age_00_04","age_05_14","age_15_64","age_65_99")
  params$n_season_cum_fit = 3 # as per RespiCompass round 1, avoid fitting the cumulative ili for an early ongoing season
  # Settings for last_year_burden
  params$last_year_burden$target = "ILIconsultationrate"
  params$last_year_burden$agegroup = "age_total"
  # Settings for arima simple
  params$arima_simple$target = "ILIconsultationrate"
  params$arima_simple$agegroup = "age_total"
  
  # ---- |-Fitting and uncertainty ----
  
  # ---- |-Flu scenarios ----
  params$scenarios =   
    tibble(
      scen_id=1:7,
      scenario_id=c("A","B","C","D","E","F","G"),
      axis_vax=c(1,1,2,2,3,3,0),
      axis_transmission=c(0,2,0,2,0,2,0)
    ) %>% left_join(
      tibble(axis_vax=c(1,2,3,0),axis_vax_name=c("opti","pess","null","status_quo")) , by = join_by(axis_vax)
    ) %>% left_join(
      tibble(axis_transmission=c(1,2,0),axis_transmission_name=c("opti","pess","status_quo")), by = join_by(axis_transmission)
    )
  
  
  
  # ---- |-Folder paths ----
  # season_cycle_round_id>-<team>-<model>.parquet (Ex. 2024_2025_1_FLU2-ISI-GLEAM.parquet)
  for (model_id in params$models_to_run) params$path_save_results[[model_id]] = paste0("./output/",
                                                                                       params$scenario_round_id,"-",
                                                                                       params$scenario_team,"-",
                                                                                       model_id,".parquet")
  
  params$path_save_figures = "./output/figures/"
  
  return(params)
}

 



