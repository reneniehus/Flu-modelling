load_flu_data = function(params=NULL){
  
  # origin of the data
  # data2 = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_ILIARIRates.csv")
  # write_csv(data2,file="./data/erviss_iliari_snapshot_2024-05-24.csv")
  # data_sentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_sentinelTestsDetectionsPositivity.csv")
  # write_csv(data_sentinel_detections,file="./data/erviss_detections_sentinel_snapshot_2024-05-24.csv")
  # data_nonsentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_nonSentinelTestsDetections.csv")
  # write_csv(data_nonsentinel_detections,file="./data/erviss_detections_nonsentinel_snapshot_2024-05-24.csv")
  data2 = read_csv(file="./data/erviss_iliari_snapshot_2024-05-24.csv")
  data2 = data2 %>% 
    mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>% 
    mutate(age = case_when(
      age == "0-4" ~ "age_00_04",
      age == "5-14" ~ "age_05_14",
      age == "15-64" ~ "age_15_64",
      age == "65+" ~ "age_65_99",
      age == "total" ~ "age_total",
      age == "unk" ~ "age_unk",
    )) %>% 
    mutate(countrycode=EU_short(countryname)) %>% 
    select( 
      country_short=countrycode,
      date=date, # see if we need to be more explicit
      target=indicator,
      agegroup=age,
      value=value
    ) 
  
  #
  data_sentinel_detections = read_csv(file="./data/erviss_detections_sentinel_snapshot_2024-05-24.csv")
  data_sentinel_detections = data_sentinel_detections %>% 
    mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>% 
    mutate(age = case_when(
      age == "0-4" ~ "age_00_04",
      age == "5-14" ~ "age_05_14",
      age == "15-64" ~ "age_15_64",
      age == "65+" ~ "age_65_99",
      age == "total" ~ "age_total",
      age == "unk" ~ "age_unk",
    )) %>% 
    mutate(country_short=EU_short(countryname)) 
  #
  data_nonsentinel_detections = read_csv(file="./data/erviss_detections_nonsentinel_snapshot_2024-05-24.csv")
  data_nonsentinel_detections = data_nonsentinel_detections %>% 
    mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>% 
    mutate(age = case_when(
      age == "0-4" ~ "age_00_04",
      age == "5-14" ~ "age_05_14",
      age == "15-64" ~ "age_15_64",
      age == "65+" ~ "age_65_99",
      age == "total" ~ "age_total",
      age == "unk" ~ "age_unk",
    )) %>% 
    mutate(country_short=EU_short(countryname)) 
  
  mout=list(
    erviss_ili_ari = data2,
    erviss_detect_sent = data_sentinel_detections,
    erviss_detect_nonsent = data_nonsentinel_detections
  )
  
  return(mout)
}


