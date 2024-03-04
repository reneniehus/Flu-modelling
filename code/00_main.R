# script that outlines the work flow for BOTH flu forecasting and scenario modelling
# ---- |-Set up ----
source("code/setup.R")

# ---- |-load task specific settings ----
source("code/settings/settings_version0.R") # changed by the user
params=settings() # calls the function that creates the params-list

# sourcing other files, models etc
source("code/load_flu_data.R")
source("code/run_flu_models.R")
source("code/process_and_save.R")
source("code/model_SIR_simple.R")
source("code/model_last_year_burden.R")

# ---- |-load flu data ----
data = load_flu_data( params ) # loads the data

# ---- |-run models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-save final output including basic figures ----
process_and_save( params, models_out ) # processing the model output, with figures and saves

# plots for checking
models_out %>% ggplot(aes(x=date,y=value,group=sample_or_quantile)) +
  geom_line(alpha=0.2) + coord_cartesian(ylim=c(0,5000))

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

source("./code/respicasting_ILI.R")
source("./code/respicasting_ARI.R")
# covid
source("./code/respicasting_covid_cases.R")
source("./code/respicasting_covid_death.R")
source("./code/respicasting_covid_hosp.R")
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





