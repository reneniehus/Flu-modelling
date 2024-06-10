# Code for model comparison

source("code/01_main_supporting/setup.R")

# ---- |-load task specific settings ----
source("code/02_settings/settings_version0.R") # changed by the user
params=settings() # calls the function that creates the params-list

# sourcing other files, models etc
source("code/01_main_supporting/load_flu_data.R")
source("code/01_main_supporting/run_flu_models.R")
source("code/01_main_supporting/process_and_save.R")
source("code/01_main_supporting/model_SIR_simple.R")
source("code/model_last_year_burden.R")

# ---- |-load flu data ----
data = load_flu_data( params ) # loads the data
data_mock = data$erviss_ili_ari %>% 
  filter(country_short == country_short_input, 
         target == params$SIR_simple$target, 
         agegroup == params$SIR_simple$agegroup) # 
data_mock %<>% filter( date%in%date_v_fit )
data_mock %>% ggplot(aes(date,value)) + geom_line()
data_mock_fit = data_mock
# Projection dates
data_mock_project = data_mock
data_mock_project$date = data_mock_project$date+365

################# Run model with age, season and vaccination status ################# 


stan_list_age2_season = list(
  n_season = 2,
  n_week_fit = nrow(data_mock_fit)*2, #
  n_day_fit = nrow(data_mock_fit)*2*7, #
  severe_obs_fit = x_obs %>% rbind(x_obs),
  season_id = rep(c(1,2), each = length(data_mock_fit$value)*7),
  n_week_project = nrow(data_mock_project),
  n_age_groups = 2,# vector
  contact_matrix = matrix(c(1,1,1,1), 2, 2),
  pop = 9e6, # real
  pop_age_group =  matrix(c(9e6,9e6)/2,2,1), # matrix
  Rnull = params$Rnull, # real
  rate_infectious = params$rate_infectious, # real
  severe_obs_notna = rep(1,nrow(data_mock_fit)*2),
  season_start = c(1,2, rep(0, nrow(data_mock_fit)*7-2), 1,2, rep(0, nrow(data_mock_fit)*7-2)),
  ve_susc = 1,
  ve_inf = 1,
  ve_severe = 1,
  delta_vax = matrix(0, 364, 2) #(x_obs %>% rbind(x_obs)) * 0
)

fit02_season_age_v =rstan::stan(
  file='./stan/SIR_multiseason_age_vax.stan',
  chains=1 ,thin=1,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list_age2_season
) 


################# Run the simple model ################# 





options(mc.cores = detectCores()-1 )
#mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script

stan_list = list(
  n_week_fit = nrow(data_mock_fit),
  severe_obs_fit = data_mock_fit$value %>% round(),
  n_week_project = nrow(data_mock_project),
  pop = 9e6,
  Rnull = params$Rnull,
  rate_infectious = params$rate_infectious
)

fit02_simple=rstan::stan(
  file='./stan/SIR_simple.stan',
  chains=8 ,thin=8,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list
) # X mins
#precis(fit02,pars = c("Rnull_eff"))

est_Rnull = precis(fit02_simple,pars = c("Rnull_eff"))
est_par = fit02_simple %>% gather_draws(SIR_ini[state],
                       prop_severe,
                       pop_infect
) %>% 
  mean_qi() -> xp; xp

################ Run the simple model ################# 


stan_list_age = list(
  n_week_fit = nrow(data_mock_fit), #
  severe_obs_fit = data_mock_fit$value %>% round %>% as.matrix(),# matrix fixme!
  n_week_project = nrow(data_mock_project),
  n_age_groups = 1,# vector
  contact_matrix =1 %>% as.matrix(),
  pop = 9e6, # real
  pop_age_group = 9e6 %>% as.matrix(), # maxtrix
  Rnull = params$Rnull, # real
  rate_infectious = params$rate_infectious # real
)

stan_list_age2 = list(
  n_week_fit = nrow(data_mock_fit), #
  severe_obs_fit = (data_mock_fit$value/2) %>% round %>% as.matrix() %>% bind_cols((data_mock_fit$value/2) %>% round %>% as.matrix()),# matrix fixme!
  n_week_project = nrow(data_mock_project),
  n_age_groups = 2,# vector
  contact_matrix = matrix(c(1,1,1,1), 2, 2),
  pop = 9e6, # real
  pop_age_group =  matrix(c(9e6,9e6)/2,2,1), # maxtrix
  Rnull = params$Rnull, # real
  rate_infectious = params$rate_infectious # real
)

x_obs = (data_mock_fit$value/2) %>% round %>% as.matrix() %>% bind_cols((data_mock_fit$value/2) %>% round %>% as.matrix()) # matrix fixme!

stan_list_age2_season = list(
  n_season = 2,
  n_week_fit = nrow(data_mock_fit)*2, #
  n_day_fit = nrow(data_mock_fit)*2*7, #
  severe_obs_fit = x_obs %>% rbind(x_obs),
  season_id = rep(c(1,2), each = length(data_mock_fit$value)*7),
  n_week_project = nrow(data_mock_project),
  n_age_groups = 2,# vector
  contact_matrix = matrix(c(1,1,1,1), 2, 2),
  pop = 9e6, # real
  pop_age_group =  matrix(c(9e6,9e6)/2,2,1), # maxtrix
  Rnull = params$Rnull, # real
  rate_infectious = params$rate_infectious, # real
  severe_obs_notna = rep(1,nrow(data_mock_fit)*2),
  season_start = c(1,2, rep(0, nrow(data_mock_fit)*7-2), 1,2, rep(0, nrow(data_mock_fit)*7-2))
)

fit02_season_age=rstan::stan(
  file='./stan/SIR_simple_multiseason_age.stan',
  chains=1 ,thin=1,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list_age2_season
) 


fit02_age=rstan::stan(
  file='./stan/SIR_age.stan',
  chains=8 ,thin=8,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list_age2
) 

est_Rnull
est_Rnull_age = precis(fit02_age,pars = c("Rnull_eff"))
est_Rnull_age
est_par
fit02_age %>% gather_draws(#SIR_init[state],
                                        prop_severe
) %>% 
  mean_qi() -> xp; xp

fit02_simple %>% gather_draws(#SIR_init[state],
  prop_severe
) %>% 
  mean_qi() -> xp; xp



######

fit02_season_age %>% spread_draws(S[t,a],I[t,a],R[t,a] ) %>% mutate(m_sum = S+I+R) %>% pull(m_sum) %>% mean()


# Run stan model for age, season and vaccination status


stan_list_age2_season = list(
  n_season = 2,
  n_week_fit = nrow(data_mock_fit)*2, #
  n_day_fit = nrow(data_mock_fit)*2*7, #
  severe_obs_fit = x_obs %>% rbind(x_obs),
  season_id = rep(c(1,2), each = length(data_mock_fit$value)*7),
  n_week_project = nrow(data_mock_project),
  n_age_groups = 2,# vector
  contact_matrix = matrix(c(1,1,1,1), 2, 2),
  pop = 9e6, # real
  pop_age_group =  matrix(c(9e6,9e6)/2,2,1), # matrix
  Rnull = params$Rnull, # real
  rate_infectious = params$rate_infectious, # real
  severe_obs_notna = rep(1,nrow(data_mock_fit)*2),
  season_start = c(1,2, rep(0, nrow(data_mock_fit)*7-2), 1,2, rep(0, nrow(data_mock_fit)*7-2)),
  ve_susc = 1,
  ve_inf = 1,
  ve_severe = 1,
  delta_vax = matrix(0, 364, 2) #(x_obs %>% rbind(x_obs)) * 0
)

fit02_season_age_v =rstan::stan(
  file='./stan/SIR_multiseason_age_vax.stan',
  chains=1 ,thin=1,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list_age2_season
) 
