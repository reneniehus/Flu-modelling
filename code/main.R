# ---- |-Set up ----
source("code/setup01.R")

# ---- |-load task specific settings ----
source("code/settings_version0.R") # changed by the user
params=settings() # calls the function that creates the params-list

# sourcing other files, models etc
source("code/load_flu_data.R")
source("code/run_flu_models.R")
source("code/process_and_save.R")

# ---- |-load flu data ----
data = load_flu_data( params ) # loads the data

# ---- |-run models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-save final output including basic figures ----
process_and_save( params, models_out ) # processing the model output, with figures and saves

