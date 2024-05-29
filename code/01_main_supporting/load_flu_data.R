load_flu_data = function( params=NULL,withforce=F, new_from_online=T ){
  
  file_doesnot_exist = !file.exists("output/data_erviss.Rdata")
  if ( file_doesnot_exist|withforce==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      erviss_ili_ari = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_ILIARIRates.csv",show_col_types = FALSE)
      data_sentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_sentinelTestsDetectionsPositivity.csv",show_col_types = FALSE)
      data_nonsentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_nonSentinelTestsDetections.csv",show_col_types = FALSE)
    }
    if (new_from_online==F) {
      # load data from local storage
      erviss_ili_ari = read_csv(file="./data/erviss_iliari_snapshot_2024-05-24.csv",show_col_types = FALSE)
      data_sentinel_detections = read_csv(file="./data/erviss_detections_sentinel_snapshot_2024-05-24.csv",show_col_types = FALSE)
      data_nonsentinel_detections = read_csv(file="./data/erviss_detections_nonsentinel_snapshot_2024-05-24.csv",show_col_types = FALSE)
    }
    
    erviss_ili_ari = erviss_ili_ari %>% 
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
    
    data_erviss=list(
      date_list_created = today(),
      erviss_ili_ari = erviss_ili_ari,
      erviss_typing_sentinel = data_sentinel_detections,
      erviss_typing_nonsentinel = data_nonsentinel_detections
    )
    save(data_erviss,file="output/data_erviss.Rdata")
  } else { load(file="output/data_erviss.Rdata") }
  
  return(data_erviss)
}


