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

# starting a forecasting script for RespiCast
source("./code/respicasting_ILI.R")
source("./code/respicasting_ARI.R")
