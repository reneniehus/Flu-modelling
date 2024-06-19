# 1: defining functions that load each data stream
# 2: define a mother function that calls each data stream function

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Defining data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load_flu_data_epi = function(data=data, params=NULL , new_from_online=T , regenerate=T ){
  
  file_doesnot_exist = !file.exists("output/epi.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      erviss_ili_ari = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_ILIARIRates.csv",show_col_types = FALSE)
      data_sentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_sentinelTestsDetectionsPositivity.csv",show_col_types = FALSE)
      data_nonsentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/snapshots/2024-05-24_nonSentinelTestsDetections.csv",show_col_types = FALSE)
    }
    if (new_from_online==F) {
      # load data from disk
      # add paths
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
    
    epi = list(
      date_list_created = today(),
      erviss_ili_ari = erviss_ili_ari,
      erviss_typing_sentinel = data_sentinel_detections,
      erviss_typing_nonsentinel = data_nonsentinel_detections
    )
    save(epi,file="output/epi.Rdata")
    
  } else { load(file="output/epi.Rdata") }
  
  # adding to data 
  data$epi = epi
  
  return(data)
}

load_flu_data_vax = function(data=data, params=NULL , new_from_online=T , regenerate=T ){
  file_doesnot_exist = !file.exists("output/vax.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
    }
    if (new_from_online==F) {
      # load data from local storage
    }
    
    vax = list(
    )
    save(vax,file="output/vax.Rdata")
    
  } else { load(file="output/vax.Rdata") }
  
  # adding to data 
  data$vax = vax
  
  return(data)
}

load_flu_data_contact = function(data=data, params=NULL , new_from_online=T , regenerate=T ){
  file_doesnot_exist = !file.exists("output/contact.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
    }
    if (new_from_online==F) {
      # load data from local storage
    }
    
    dat_contact = list(
    )
    save(dat_contact,file="output/contact.Rdata")
    
  } else { load(file="output/contact.Rdata") }
  
  # adding to data 
  data$contact = dat_contact
  
  return(data)
}

load_flu_data_demography = function(data=data, params=NULL , new_from_online=F , regenerate=T ){
  file_doesnot_exist = !file.exists("output/demography.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      # required libraries
      source('db/logger.R')
      source('db/sql_utils.R')
      # Logger is needed for running SQL utils
      logger <- forge_logger()(logLevel = 'INFO')
      # This is where the SQL Profiles are stored (always needed)
      dbDir <- 'db/'
      # This is where the SQL templates are stored (needed only if you're using them)
      SQLTemplatePath <- 'db/templates/sql/'
      ## Querying simple data ----
      pop_data <- read_data(table = 'out.DM_Population_ByCountryEU',
                            connInfo = 'pop')
      pop_data %>% as_tibble() %>% 
        filter(ReportYear==2024,CountryCode%in%countries_short) %>% 
        filter(AgeGroup %in% c("Age00_04",
                               "Age05_09",
                               "Age10_14",
                               "Age15_19", # up to 16 years, single age brackets exist, too
                               "Age20_24",
                               "Age25_29",
                               "Age30_34",
                               "Age35_39",
                               "Age40_44",
                               "Age45_49",
                               "Age50_54",
                               "Age55_59",
                               "Age60_64",
                               "Age65_69",
                               "Age70_74",
                               "Age75_79",
                               "Age80_84",
                               "Age85_89",
                               "Age90_94",
                               "Age95+") ) %>% 
        group_by(country=CountryCode,age_group=AgeGroup) %>% 
        mutate(Population=as.numeric(Population)) %>% 
        summarise(population=sum(Population)) %>% ungroup() -> mdat
      write_fst(mdat,path="output/population_pyramid.fst")
    }
    if (new_from_online==F) {
      # load data from local storage
      mdat = read_fst(path="output/population_pyramid.fst") %>% as_tibble()
    }
    
    dat_demography = list(
      population_pyramid = mdat
    )
    save(dat_demography,file="output/demography.Rdata")
    
  } else { load(file="output/demography.Rdata") }
  
  # adding to data 
  data$demography = dat_demography
  
  return(data)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Mother function: calling the data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load_flu_data = function( params=NULL , new_from_online=T, regenerate=F ){
  
  data = list() # reset data list
  
  data = load_flu_data_epi( data=data, params=NULL , new_from_online=T , regenerate=regenerate)
  
  data = load_flu_data_vax( data=data, params=NULL , new_from_online=F , regenerate=regenerate)
  
  data = load_flu_data_contact( data=data, params=NULL , new_from_online=F , regenerate=regenerate)
  
  data = load_flu_data_demography( data=data, params=NULL , new_from_online=F , regenerate=regenerate)
  
  return(data)
}

