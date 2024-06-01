# script that outlines the work flow for BOTH flu forecasting and scenario modelling
# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R") # changed by the user
params=settings() # calls the function that creates the params-list

# sourcing other files, models etc
source("code/01_main_supporting/flu_functions.R")
source("code/01_main_supporting/load_flu_data.R")
source("code/01_main_supporting/run_flu_models.R")
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/model_SIR_simple.R")
source("code/01_main_supporting/model_last_year_burden.R")

# ---- |-load flu data ----
data = load_flu_data( params ) # loads the data

# ---- |-run models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-save final output including basic figures ----
process_and_save( params, models_out ) # processing the model output, with figures and saves

# ---- |- Run special analyses
source("code/03_special_analyses/season_country_var_rt.R")
source("code/03_special_analyses/erviss_data_look.R")
source("code/03_special_analyses/exploring_SIR/SIR_paras_explore.R") # explore what SIR paras do
if (F) source("code/03_special_analyses/forecasting/norrsken.R") # forecasting modelling

# temporary code for quick checking
# ask Leah/Nick about overlap of sent and nonsent typing data [Asked Leah]
# fit ALL 25 countries: ILI
# fit ALL 25 countries: ILI*sent_typing
# eyeballing method: identify good and poor data/fits
# unless all fits look amazing: fit ALL 25 countries: ILI*nonsent_typing
# -> compare fitted parameters  ILI*sent_typing versus ILI*nonsent_typing

# fit all poor-fit countries: ILI*nonsent_typic (or combination of sent&nonsent)
models_out$multiseason$p
