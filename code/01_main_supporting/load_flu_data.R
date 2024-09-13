# 1: defining functions that load each data stream
# 2: define a mother function that calls each data stream function

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Defining data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load_flu_data_epi = function(data=data, params=NULL, regenerate=T , new_from_online=T){
  
  file_doesnot_exist = !file.exists("output/epi.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      pr=paste("Loading epi data from github ... \n"); cat(green(pr))
      erviss_ili_ari = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/ILIARIRates.csv",show_col_types = FALSE)
      data_sentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/sentinelTestsDetectionsPositivity.csv",show_col_types = FALSE)
      data_nonsentinel_detections = read_csv(file="https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/nonSentinelTestsDetections.csv",show_col_types = FALSE)
      data_respicompass_iliplus = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/target-data/influenza/ili_plus.csv",show_col_types = FALSE)
      # save locally
      erviss_ili_ari %>% write_csv(file="./data/erviss_iliari_snapshot_2024-05-24.csv")
      data_sentinel_detections %>% write_csv(file="./data/erviss_detections_sentinel_snapshot_2024-05-24.csv")
      data_nonsentinel_detections %>% write_csv(file="./data/erviss_detections_nonsentinel_snapshot_2024-05-24.csv")
      data_respicompass_iliplus %>% write_csv(file="./data/data_respicompass_iliplus.csv")
    }
    if (new_from_online==F) {
      pr=paste("Loading epi data from disk ... \n"); cat(green(pr))
      # load data from local storage
      erviss_ili_ari = read_csv(file="./data/erviss_iliari_snapshot_2024-05-24.csv",show_col_types = FALSE)
      data_sentinel_detections = read_csv(file="./data/erviss_detections_sentinel_snapshot_2024-05-24.csv",show_col_types = FALSE)
      data_nonsentinel_detections = read_csv(file="./data/erviss_detections_nonsentinel_snapshot_2024-05-24.csv",show_col_types = FALSE)
      data_respicompass_iliplus = read_csv(file="./data/data_respicompass_iliplus.csv",show_col_types = FALSE)
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
    
    #
    respicompass_iliplus = data_respicompass_iliplus %>% mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>% 
      mutate(age = case_when(
        age == "0-4" ~ "age_00_04",
        age == "5-14" ~ "age_05_14",
        age == "15-64" ~ "age_15_64",
        age == "65+" ~ "age_65_99",
        age == "total" ~ "age_total",
        age == "unk" ~ "age_unk",
      )) %>% 
      mutate(countrycode=EU_short(location_name),target="ili_plus") %>% 
      select( 
        country_short=countrycode,
        date=date, # see if we need to be more explicit
        target,
        agegroup=age,
        value=value
      ) 
    epi = list(
      date_epilist_created = today(),
      erviss_ili_ari = erviss_ili_ari,
      erviss_typing_sentinel = data_sentinel_detections,
      erviss_typing_nonsentinel = data_nonsentinel_detections,
      respicompass_iliplus = respicompass_iliplus
    )
    save(epi,file="output/epi.Rdata")
    
  } else { load(file="output/epi.Rdata") }
  
  # adding to data 
  data$epi = epi
  
  return(data)
}

load_flu_data_vax = function(data=data, params=NULL , regenerate=T  , new_from_online=T){
  file_doesnot_exist = !file.exists("output/vax.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      data_vax = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/influenza/vaccination/influenza_vax_scenarios.csv",show_col_types = FALSE)
      data_vax_hist = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/influenza/vaccination/vaccine_coverage_65plus.csv",show_col_types = FALSE)
      data_vax_hist_all = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/influenza/vaccination/vaccine_coverage_all.csv",show_col_types = F)
      
      # write
      data_vax %>% write_csv(file="data/vax_flu_scenarios.csv")
      data_vax_hist %>% write_csv(file="data/vax_flu_history.csv")
      data_vax_hist_all %>% write_csv(file="data/vax_flu_history_all.csv")
      
    }
    if (new_from_online==F) {
      # load data from local storage
      data_vax = read_csv(file="data/vax_flu_scenarios.csv",show_col_types = F )
      data_vax_hist = read_csv(file="data/vax_flu_history.csv",show_col_types = F )
      data_vax_hist_all = read_csv(file="data/vax_flu_history_all.csv",show_col_types = F )
      
    }
    
    vax = list(
      data_vax = data_vax %>% mutate(vaccine_coverage=vaccine_coverage/100) %>% pivot_wider(names_from = "scenario", values_from = vaccine_coverage),
      data_vax_history = data_vax_hist %>% mutate(vaccine_coverage=as.numeric(vaccine_coverage)/100 ) %>% mutate(season = str_replace(season, "-", "/")),
      data_vax_history_all = data_vax_hist %>% mutate(vaccine_coverage=as.numeric(vaccine_coverage)/100 ) %>% mutate(season = str_replace(season, "-", "/"))
    )
    save(vax,file="output/vax.Rdata")
    
  } else { load(file="output/vax.Rdata") }
  
  # adding to data 
  data$vax = vax
  
  return(data)
}

load_flu_data_contact = function(data=data, params=NULL , regenerate=T  , new_from_online=F){
  file_doesnot_exist = !file.exists("output/contact.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      stop("Data available only locally")
      
    }
    if (new_from_online==F) {
      # load data from local storage
      xdata = list()
      xlocations = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/supporting-files/locations_iso2_codes.csv",show_col_types = F)
      for (country_i in xlocations$location_name){ # country_i = xlocations$location_name[1]
        contacts = 0
        try({ # As there are two files, each with half the countries, try to get contact matrix for a given country from both files
          contacts = read_excel("./data/MUestimates_all_locations_1.xlsx", sheet = country_i, col_names = F, .name_repair = "unique_quiet", skip = 1)
        }, silent = T)
        try({
          contacts = read_excel("./data/MUestimates_all_locations_2.xlsx", sheet = country_i, col_names = F, .name_repair = "unique_quiet")
        }, silent = T)
        xdata[[country_i]] = contacts
        
      }
    }
    
    dat_contact = xdata
    save(dat_contact,file="output/contact.Rdata")
    
  } else { load(file="output/contact.Rdata") }
  
  # adding to data 
  data$contact = dat_contact
  
  return(data)
}

load_flu_data_helpers_respicompass = function(data=data, params=NULL , regenerate=T  , new_from_online=T){
  
  file_doesnot_exist = !file.exists("output/respicompass_helpers.Rdata")
  
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      xlocations = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/supporting-files/locations_iso2_codes.csv",show_col_types = F)
      xweeks = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/supporting-files/iso_weeks.csv",show_col_types = F)
      # write to disk
      xlocations %>% write_csv(file="output/respicompass_locations.csv")
      xweeks %>% write_csv(file="output/respicompass_weeks.csv")
    }
    if (new_from_online==F) {
      # load data from local storage
      xlocations = read_csv(file="output/respicompass_locations.csv",show_col_types = F)
      xweeks = read_csv(file="output/respicompass_weeks.csv",show_col_types = F)
    }
    
    dat_helpers = list(
      iso2_code = xlocations,
      iso_weeks = xweeks
    )
    save(dat_helpers,file="output/respicompass_helpers.Rdata")
    
  } else { load(file="output/respicompass_helpers.Rdata") }
  
  # adding to data 
  data$helpers_respicompass = dat_helpers
  
  return(data)
}

load_flu_data_demography_ECDC = function(data=data, params=NULL , regenerate=T  , new_from_online=T){
  file_doesnot_exist = !file.exists("output/demography.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      # required libraries
      pr=paste("Loading demography data from database ... \n"); cat(green(pr))
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
      write_fst(mdat,path="data/population_pyramid.fst")
    }
    if (new_from_online==F) {
      # load data from local storage
      pr=paste("Loading demography data from disk ... \n"); cat(green(pr))
      mdat = read_fst(path="data/population_pyramid.fst") %>% as_tibble()
    }
    
    dat_demography = list(
      population_pyramid = mdat
    )
    save(dat_demography,file="output/demography.Rdata")
    
  } else { load(file="output/demography.Rdata") }
  
  # adding to data 
  data$demography_ECDC = dat_demography
  
  return(data)
}

load_flu_data_demography_respicast = function(data=data, params=NULL , regenerate=T  , new_from_online=T){
  file_doesnot_exist = !file.exists("output/demography_respicast.Rdata")
  if ( file_doesnot_exist|regenerate==T ) {
    
    if (new_from_online==T) {
      # load data freshly from the internet
      pop_df = NULL
      pop_fine_df = NULL
      xlocations = read_csv(file="output/respicompass_locations.csv",show_col_types = F)
      country_v = xlocations$location_name
      for (country_i in country_v) {
        pr=paste("> Loading pop data for:",country_i,"... \n"); cat(green(pr))
        read_file=paste0("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/miscellaneous/population/",country_i,"_aggr.csv")
        read_file_fine = paste0("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/miscellaneous/population/",country_i,".csv")
        xdf = read_csv(read_file,show_col_types = FALSE)
        xdf_fine = read_csv(read_file_fine,show_col_types = FALSE)
        xdf$country = country_i
        xdf_fine$country = country_i
        pop_df=rbind(pop_df,xdf)
        pop_fine_df=rbind(pop_fine_df,xdf_fine)
      }
      pop_df = pop_df %>% select(country,age_group,population)
      pop_df %>% write_csv("output/population_pyramid_respicast.csv")
      
      pop_fine_df = pop_fine_df %>% select(country,age_group,population)
      pop_fine_df %>% write_csv("output/population_pyramid_fine_respicast.csv")
      
    }
    if (new_from_online==F) {
      # load data from local storage
      pr=paste("Loading respicast demography data from disk ... \n"); cat(green(pr))
      pop_df = read_csv("output/population_pyramid_respicast.csv",show_col_types = F)
      pop_fine_df = read_csv("output/population_pyramid_fine_respicast.csv",show_col_types = F)
    }
    
    dat_demography = list(
      population_pyramid = pop_df,
      population_pyramid_fine = pop_fine_df
    )
    save(dat_demography,file="output/demography_respicast.Rdata")
    
  } else { load(file="output/demography_respicast.Rdata") }
  
  # adding to data 
  data$demography_respicast = dat_demography
  
  return(data)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Mother function: calling the data-loading functions for each data stream ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load_flu_data = function( params=NULL , regenerate=F,new_from_online=T  ){
  
  data = list() # reset data list
  
  data = load_flu_data_epi( data=data, params=NULL , new_from_online=T , regenerate=T )
  
  data = load_flu_data_vax( data=data, params=NULL , new_from_online=new_from_online , regenerate=regenerate)
  
  data = load_flu_data_contact( data=data, params=NULL , new_from_online=F , regenerate=regenerate)
  
  data = load_flu_data_helpers_respicompass( data=data, params=NULL , new_from_online=new_from_online , regenerate=regenerate)
  
  data = load_flu_data_demography_ECDC( data=data, params=NULL , new_from_online=F , regenerate=F)
  
  data = load_flu_data_demography_respicast( data=data, params=NULL , new_from_online=new_from_online , regenerate=regenerate)
  
  return(data)
}

