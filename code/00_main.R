# ---- |-Start empty ----
rm(list = ls(all.names = TRUE)); gc() # clear environment & memory

# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R"); params=settings() # settings_version_X.R script to be changed by high-level user

# ---- |-sourcing support scripts ----
source("code/01_main_supporting/flu_functions.R")
source("code/01_main_supporting/load_flu_data.R")
source("code/01_main_supporting/run_flu_models.R")
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/send_report.R")
source("code/01_main_supporting/model_SIR_multiseason.R")
source("code/06_sandbox/generate_ili_epi_test.R")

# ---- |-load flu data ----
data = load_flu_data( params, regenerate = F, new_from_online = F) # loads the data # regenerate=T recreates the data lists, new_from_online=T uses the online versions for recreation

# ---- |-run flu models (fitting and projections) ----
models_out = run_flu_models( params, data ) # runs the model scripts
# ---- |-process and save model output ----
rep_list = process_and_save( params, data, models_out, save_submission=params$save_submission );  # processing the model output, with figures and saves

# ---- |-report ----
if (F) rmarkdown::render("code/03_report/report_overview.Rmd") # requires: params, data , rep_list

# ---- |-Send ----
if (F) send_report(params)

# ---- |- Run special analyses
if (F) source("code/04_special_analyses/burden_vs_vax.R") # can be run after data is loaded (load_flu_data)
if (F) source("code/04_special_analyses/ili_burden_over_seasons.R") # can be run after data is loaded (load_flu_data)
if (F) source("code/04_special_analyses/ili_vs_sari_burden.R") # can be run after data is loaded (load_flu_data)
if (F) source("code/04_special_analyses/season_country_var_rt.R")
if (F) source("code/04_special_analyses/erviss_data_look.R")
if (F) source("code/04_special_analyses/exploring_SIR/SIR_paras_explore.R") # explore what SIR paras do
if (F) source("code/04_special_analyses/forecasting/norrsken.R") # forecasting modelling

# ---- |- The end

# (temporary code for any quick checking)