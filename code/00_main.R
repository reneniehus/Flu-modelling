# main script that outlines the work flow for BOTH flu forecasting and scenario modelling

# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R") # this script is what is changed by a high-level user
params=settings() # creates the params-list
figs = list() # create an empty list for figures

# sourcing other files, models etc
source("code/01_main_supporting/flu_functions.R")
source("code/01_main_supporting/load_flu_data.R")
source("code/01_main_supporting/run_flu_models.R")
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/model_SIR_simple.R")
source("code/01_main_supporting/model_last_year_burden.R")

# ---- |-load flu data ----
data = load_flu_data( params, regenerate = F, new_from_online = F) # loads the data
# regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation

# ---- |-run models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-save final output including basic figures ----
process_and_save( params, models_out ) # processing the model output, with figures and saves

# ---- |- Run special analyses
if (F) source("code/03_special_analyses/season_country_var_rt.R")
if (F) source("code/03_special_analyses/erviss_data_look.R")
if (F) source("code/03_special_analyses/exploring_SIR/SIR_paras_explore.R") # explore what SIR paras do
if (F) source("code/03_special_analyses/forecasting/norrsken.R") # forecasting modelling

# (temporary code for any quick checking)


