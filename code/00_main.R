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
# bad fit: AT, IT


if (F){
  df = NULL
  for (i in 1:length(models_out$other)) {
    x = models_out$other[[i]]$pars_df
    x = as_tibble(x, rownames = "para")
    x$country = names(models_out$other)[i]
    x$country_long = names(models_out$other)[i] %>% EU_long()
    df = rbind(df,x)
  }
  df %>% filter(para=="prop_ili_mu") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) + 
    scale_x_log10()
  df %>% filter(para=="prop_ili_mu") %>% pull(mean) %>% 
    
    df %>% filter(para=="SIR_ini_mu[1]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`))
  
  df %>% filter(para=="SIR_ini_mu[2]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10()
}

