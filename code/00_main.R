rm(list = ls(all.names = TRUE)) # clear environment
gc() # clear memory 
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
source("code/01_main_supporting/model_SIR_multiseason.R")
source("code/04_sandbox/generate_ili_epi_test.R")

# ---- |-load flu data ----
data = load_flu_data( params, regenerate = F, new_from_online = F) # loads the data
# regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation

# ---- |-run flu models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-save output ----
if (F) process_and_save( params, models_out ) # processing the model output, with figures and saves

# ---- |- Run special analyses
if (F) source("code/03_special_analyses/season_country_var_rt.R")
if (F) source("code/03_special_analyses/erviss_data_look.R")
if (F) source("code/03_special_analyses/exploring_SIR/SIR_paras_explore.R") # explore what SIR paras do
if (F) source("code/03_special_analyses/forecasting/norrsken.R") # forecasting modelling

# (temporary code for any quick checking)
# models_out$figs_prefit$fit_seasons_countries
# models_out$other$IT$plot_fit
models_out$other$DK$plot_fit


# 4 age groups
stan_list$contact_matrix = matrix(data=rep(4,each=16),nrow=4)
stan_list = generate_ili_epi_test(par = c(NA),stan_list)
stan_list$ili_obs_fit %>% sum() # 6376188 (2 age groups, no noise)

# 2 age groups
stan_list = make_stan_list(params,data,all_season_fit_wide,country_short_input,vax_country,pop_country,age_collapse="two")
stan_list$contact_matrix = matrix(data=c(1/2,1/2,1/2,1/2),nrow=2)
stan_list$delta_vax = stan_list$delta_vax*0 
stan_list$pop_age_group = stan_list$pop_age_group*0 + stan_list$pop/2
stan_list = generate_ili_epi_test(par = c(NA),stan_list)
stan_list$n_age_groups
stan_list$contact_matrix
x2 = stan_list$ili_obs_fit
stan_list$ili_obs_fit %>% sum() # 8385366

# 1 age groups
stan_list = make_stan_list(params,data,all_season_fit_wide,country_short_input,vax_country,pop_country,age_collapse="one")
stan_list$delta_vax = stan_list$delta_vax*0 
stan_list$pop_age_group
stan_list = generate_ili_epi_test(par = c(NA),stan_list)
stan_list$n_age_groups
stan_list$contact_matrix
x1 = stan_list$ili_obs_fit
stan_list$ili_obs_fit %>% sum() # 8385360

