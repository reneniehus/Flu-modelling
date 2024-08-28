

source("code/01_main_supporting/setup.R")

# ---- |-Respicasting: Norrsken ----
# starting a forecasting script for RespiCast
# load pop data
#axiliary data

# ERVISS: week 52 now there (8 Jan 2024)
# COVID forecastdata updated (https://github.com/european-modelling-hubs/covid19-forecast-hub-europe/tree/main/data-truth/ECDC)
# flu data updated 
# ari data updated

aux = list()
aux$pops = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-locations/locations_eu.csv",
                    show_col_types=F)

source("./code/04_special_analyses/forecasting/respicasting_ILI.R")
source("./code/04_special_analyses/forecasting/respicasting_ARI.R")
# covid
source("./code/04_special_analyses/forecasting/respicasting_covid_cases.R")
source("./code/04_special_analyses/forecasting/respicasting_covid_death.R")
source("./code/04_special_analyses/forecasting/respicasting_covid_hosp.R")
# merge: case,death,hosp
if (T){
  # forecast_date
  # target
  # target_end_date
  # location
  # type
  # quantile
  # value
  
  # log / green
  mod_cases = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_greencase.csv"),col_types = cols(.default = "c"))
  mod_death = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_greendeath.csv"),col_types = cols(.default = "c"))
  mod_hosp = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_greenhosp.csv"),col_types = cols(.default = "c"))
  
  #mod_death = mod_death %>% filter(!(target=="4 wk ahead inc death"))
  
  
  mod_combined = bind_rows(mod_cases,mod_death,mod_hosp)
  # inspect by hand
  g(mod_combined)
  clc()
  describe(mod_combined)
  # save 
  write_csv( mod_combined , file=paste0("./output/covid-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_green.csv") )
  
  
  # sqrt/ blue
  mod_cases = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_bluecase.csv"),col_types = cols(.default = "c"))
  mod_death = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_bluedeath.csv"),col_types = cols(.default = "c"))
  mod_hosp = read_csv(file=paste0("./output/covid-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_bluehosp.csv"),col_types = cols(.default = "c"))
  
  mod_combined = bind_rows(mod_cases,mod_death)
  # inspect by hand
  g(mod_combined)
  clc()
  describe(mod_combined)
  
  # save 
  write_csv( mod_combined , file=paste0("./output/covid-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_blue.csv") )
  
}
