source("./code/setup01.R")
source("./code/Flu/DK_model/bayesian_functions.R")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Data ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Flu data from Australia: https://www.health.gov.au/resources/publications/nndss-public-dataset-influenza-laboratory-confirmed
# df_aust = read_xlsx(file="../COVID19_modelling/data/national-notifiable-diseases-surveillance-system-nndss-public-dataset-influenza-laboratory-confirmed-dataset.xlsx")
# instead
# https://www.who.int/tools/flunet
# df_aust = read_xlsx(path="../COVID19_modelling/data/DataExport_101222_2.xlsx")
# #df_aust = read_csv("../COVID19_modelling/data/data_AUS_ILI.csv")
# df_aust$ISO_YEAR %>% range()
# 
# 
# 
# df_aust %>% 
#   filter(`COUNTRY/AREA/TERRITORY`=="Denmark",ORIGIN_SOURCE=="SENTINEL") %>% 
#   ggplot(aes(ISO_SDATE,as.numeric(INF_ALL))) + geom_point()
# df_aust %>% 
#   filter(`WHOREGION`=="EUR",ORIGIN_SOURCE=="SENTINEL") %>% 
#   ggplot(aes(ISO_SDATE,as.numeric(INF_ALL))) + geom_point() + 
#   facet_wrap(~`COUNTRY/AREA/TERRITORY`)
# 
# df_aust %>% 
#   filter(`COUNTRY/AREA/TERRITORY`=="Denmark",ORIGIN_SOURCE=="SENTINEL") %>% 
#   select(ISO_SDATE,INF_ALL) %>% arrange(ISO_SDATE) -> df_stan
# pop_denmark = 5.8 * 1e6

# data for Australia from Rok
# data_AUS_ILI
aus_ili = read_csv(file="./data/data_AUS_ILI.csv")
aus_ili %>% filter(date>"2022-03-01") %>% ggplot(aes(x=date,y=value)) + geom_line()
aus_ili %>% filter(date>"2022-03-01") -> df_aus

#aus_ili %>% filter(date>"2016-03-01", date<"2016-12-01") %>% ggplot(aes(x=date,y=value)) + geom_line()
#aus_ili %>% filter(date>"2016-03-01", date<"2016-12-01") -> df_aus


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Priors ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# some basic prior explorations
source("./code/Rene sandbox/priors_explore.R")
# current prior exploration

# prob severe, but consider that infections here are "fully immunising" infections
rnorm( 3000,mean=logit(0.02), sd=abs(logit(0.02)-logit(0.20)) ) %>% inv_logit() %>% dens(show.HPDI = 0.8)

# S ini for flu # genetic drift/shift of virus + waning immunity since last waves/vaccination
rnorm( 3000,mean=logit(0.7), sd=abs(logit(0.7)-logit(0.4)) ) %>% inv_logit() %>% dens(show.HPDI = 0.8)
# I _ini
rnorm( 3000,mean=-10, sd=abs((-10)-(-15) ) ) %>% dens(show.HPDI = 0.8)
inv_logit(-10)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###  Main projects ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Plots, Poteriors, checks ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Speciality scripts ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("./code/flu/DK_model/SIR_paras_explore.R")
source("./code/flu/DK_model/fit_project_flu_sims.R")
source("./code/flu/DK_model/fit_project_flu_real.R")
source("./code/flu/DK_model/fit_diff_length.R")

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Modelling ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# set stan options
# options(mc.cores = detectCores()-1 )
# rstan_options(auto_write = TRUE)
# parallel:::setDefaultClusterOptions(setup_strategy = "parallel") # not sequential
# # fit 1
# stan_list = list(
#   n = nrow(df_hosp),
#   hosp_obs_log = log(df_hosp$hosp_obs),
#   pop=df_hosp$pop[1]
# )
# myinitofit <- lapply(1:4, function(id) list( beta_logit=-0.847,
#                                              prop_hosp_logit=-3.891,
#                                              S_ini_logit=1.38,
#                                              simga=0.02) )
# 
# start_time <- Sys.time()
# fit <- stan(file='./stan/immu_dyn_01_vocs.stan',data=stan_list,seed=14,
#             iter=800,init=myinitofit,chains=4)
# end_time <- Sys.time()
# end_time - start_time -> fitting_duration # duration for entire fit
# 
# # super quick check
# precis(fit,depth = 1) #,pars="eff_RM_EU"
# precis(fit_B,depth = 2)
# summary(fit, c('a', 'phi'))$summary
# #
# get_variables(fit)
# 
# fit %>%
#   gather_draws(mu_bx1,mu_bx2) %>%
#   mean_qi()
# fit_B %>%
#   gather_draws(g) %>%
#   mean_qi()
# 
# ggplot(datplot, aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
#   geom_pointrange() +
#   geom_vline(xintercept = 0, linetype = 'dotted')