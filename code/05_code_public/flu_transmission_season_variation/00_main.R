# ---- |-Set working directory ----
# user to set the path where the 00_main.R script is stored
path_of_00_main_file = "~/R modelling/Flu-modelling/code/05_code_public/season_variance_transmission"
# working directory is set to the folder containing 00_main.R 
setwd(path_of_00_main_file)

# ---- |-Set up ----
source("main_supporting/setup.R")

# ---- |-load task specific settings ----
source("main_supporting/settings_version0.R") # this script is what is changed by a high-level user
params=settings() # creates the params-list

# sourcing required functions, models etc
source("main_supporting/flu_functions.R")
source("main_supporting/load_flu_data.R")
source("main_supporting/run_flu_models.R")
source("main_supporting/model_SIR_simple.R")

# ---- |-load flu data ----
data = load_flu_data( params )

# ---- |-run models or analyses ----
models_out = run_flu_models( params, data ) 

# ---- |-run analyses on model results ----
source( "main_supporting/season_country_var_rt.R" ) # processing the model output, with figures and saves



