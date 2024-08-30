rm(list = ls(all.names = TRUE)) # clear environment
gc() # clear memory 
t1 <- Sys.time()

# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R") # this script is what is changed by a high-level user
params=settings() # creates the params-list

# sourcing other files, models etc
source("code/01_main_supporting/flu_functions.R")
source("code/01_main_supporting/load_flu_data.R")
source("code/01_main_supporting/run_flu_models.R")
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/send_report.R")
source("code/01_main_supporting/model_SIR_multiseason.R")
source("code/06_sandbox/generate_ili_epi_test.R")

# ---- |-load flu data ----
data = load_flu_data( params, regenerate = F, new_from_online = F) # loads the data
# regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation

# ---- |-run flu models (i.e. fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts

# ---- |-process and save model output ----
rep_list = process_and_save( params, data, models_out, save_submission=params$save_submission ) # processing the model output, with figures and saves

# ---- |-report ----
rmarkdown::render("code/03_report/report_overview.Rmd") # requires: params, data , rep_list

# ---- |-Send ----
send_report()

# ---- |- Run special analyses
if (F) source("code/04_special_analyses/season_country_var_rt.R")
if (F) source("code/04_special_analyses/erviss_data_look.R")
if (F) source("code/04_special_analyses/exploring_SIR/SIR_paras_explore.R") # explore what SIR paras do
if (F) source("code/04_special_analyses/forecasting/norrsken.R") # forecasting modelling
t2 <- Sys.time(); t2-t1 # 20 mins

# (temporary code for any quick checking)


