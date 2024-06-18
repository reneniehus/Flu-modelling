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
data_mock = data$epi$erviss_ili_ari %>% 
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


######### Contact matrix and demography data

library("socialmixr")
library("eurostat")

data(polymod)
contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0, 1, 5, 15))
survey_countries(polymod)
c <- contact_matrix(survey = polymod, countries = "Germany", age.limits = c(0, 60), symmetric = TRUE, per.capita = TRUE)
c$matrix.per.capita %*% diag(c$demography$population) 

c$matrix.per.capita %*% diag() 



obtain_demography = function( country, ten_year_brackets = FALSE ){
  path_core_functions <- "./"
  if ( ten_year_brackets ){
    if ( country == "Australia" ){
      #Source:  
      #https://www.abs.gov.au/statistics/people/population/regional-population-age-and-sex/2021#data-download
      #For 5 year age brackets, 0-4yr, 5-9yr, ..., 80-84yr, 85+yr
      agegroup_sizes <- c(1509959, 
                          1616774,
                          1623892,
                          1479632,	
                          1623384,	
                          1822031,	
                          1899620,	
                          1867387,	
                          1654500,	
                          1650035,	
                          1611554,	
                          1550507,	
                          1465025,	
                          1280143,
                          1146773,	
                          807195,	
                          545408,	
                          534260 )
      agegroup_sizes <- c( "0-09yr"= sum( agegroup_sizes[ 1:2 ] ),
                           "10-19yr"= sum( agegroup_sizes[ 3:4 ] ),
                           "20-29yr"= sum( agegroup_sizes[ 5:6 ] ),
                           "30-39yr"= sum( agegroup_sizes[ 7:8 ] ),
                           "40-49yr"= sum( agegroup_sizes[ 9:10 ] ),
                           "50-59yr"= sum( agegroup_sizes[ 11:12 ] ),
                           "60-69yr"= sum( agegroup_sizes[ 13:14 ] ),
                           "70-79yr"= sum( agegroup_sizes[ 15:16 ] ),
                           "80+yr"= sum( agegroup_sizes[ 17:18 ] ) )
    } else {
      
      suppressMessages(
        df_demography <- read_excel( paste0( path_core_functions, "data/eurostat_demography_data_5year_brackets.xlsx" ), sheet=3 )
      )
      colnames( df_demography )[1] <- "A"
      df_demography %<>% filter( grepl( "Prop", A ) | grepl( "GEO", A ) )
      df_demography <- df_demography[ , which( !is.na( df_demography[1,] ) ) ]
      colnames( df_demography )[1] <- "agegroup"
      colnames( df_demography )[2:ncol(df_demography)] <- df_demography[1,2:ncol(df_demography)]
      df_demography <- df_demography[2:nrow(df_demography),]
      colnames( df_demography )[grepl( "Kosovo", colnames( df_demography ) )] <- "Kosovo"
      colnames( df_demography )[grepl( "Germany", colnames( df_demography ) )] <- "Germany"
      colnames( df_demography )[grepl( "European Union", colnames( df_demography ) )] <- "EU"
      #5-9 years row is in the wrong place, reorder rows:
      df_demography <- bind_rows( df_demography[1,], df_demography[12,], df_demography[2:11,], df_demography[13:17,] )
      df_demography$agegroup <- str_remove(df_demography$agegroup, "Proportion of population aged ")
      df_demography$agegroup <- str_replace(df_demography$agegroup, " years", "yr")
      df_demography$agegroup[ df_demography$agegroup=="80yr and more"] <- "80+yr"
      agegroup_distribution <- as.numeric( df_demography[[country]] )
      agegroup_distribution <- c( "0-09yr"=sum( agegroup_distribution[ 1:2 ] ),
                                  "10-19yr"=sum( agegroup_distribution[ 3:4 ] ),
                                  "20-29yr"=sum( agegroup_distribution[ 5:6 ] ),
                                  "30-39yr"=sum( agegroup_distribution[ 7:8 ] ),
                                  "40-49yr"=sum( agegroup_distribution[ 9:10 ] ),
                                  "50-59yr"=sum( agegroup_distribution[ 11:12 ] ),
                                  "60-69yr"=sum( agegroup_distribution[ 13:14 ] ),
                                  "70-79yr"=sum( agegroup_distribution[ 15:16 ] ),
                                  "80+yr"= agegroup_distribution[ 17 ] )
      agegroup_distribution <- agegroup_distribution/sum(agegroup_distribution)
    }
  }else{
    #df_demography <- read_excel( paste0( path_core_functions, "data/eurostat_demography_data.xlsx" ), sheet=3 )
    df_demography <- read_fst( paste0( path_core_functions, "data/eurostat_demography_data.fst" ) )
    colnames( df_demography )[1] <- "A"
    df_demography %<>% filter( grepl( "Prop", A ) | grepl( "GEO", A ) )
    df_demography <- df_demography[ , which( !is.na( df_demography[1,] ) ) ]
    colnames( df_demography )[1] <- "agegroup"
    colnames( df_demography )[2:ncol(df_demography)] <- df_demography[1,2:ncol(df_demography)]
    df_demography <- df_demography[2:nrow(df_demography),]
    colnames( df_demography )[grepl( "Kosovo", colnames( df_demography ) )] <- "Kosovo"
    colnames( df_demography )[grepl( "Germany", colnames( df_demography ) )] <- "Germany"
    colnames( df_demography )[grepl( "European Union", colnames( df_demography ) )] <- "EU"
    #5-9 years row is in the wrong place, reorder rows:
    df_demography <- bind_rows( df_demography[1,], df_demography[8,], df_demography[2:7,], df_demography[9:13,] )
    
    df_demography$agegroup <- str_remove(df_demography$agegroup, "Proportion of population aged ")
    df_demography$agegroup <- str_replace(df_demography$agegroup, " years", "yr")
    df_demography$agegroup[ df_demography$agegroup=="80yr and more"] <- "80+yr"
    agegroup_distribution <- as.numeric( df_demography[[country]] )
    agegroup_distribution <- c( "0-04yr"=agegroup_distribution[ 1 ],
                                "05-09yr"=agegroup_distribution[ 2 ],
                                "10-14yr"=agegroup_distribution[ 3 ],
                                "15-17yr"=agegroup_distribution[ 4 ]*3/5, #"15-17yr" is approximately 3/5*"15-19yr"
                                "18-24yr"=agegroup_distribution[ 4 ]*2/5 + agegroup_distribution[ 5 ], #"19-24yr" is approximately 2/5*"15-19yr"+"20-24yr" 
                                "25-49yr"=agegroup_distribution[ 6 ],
                                "50-59yr"=agegroup_distribution[ 7 ] + agegroup_distribution[ 8 ],
                                "60-69yr"=agegroup_distribution[ 9 ] + agegroup_distribution[ 10 ],
                                "70-79yr"=agegroup_distribution[ 11 ] + agegroup_distribution[ 12 ],
                                "80+yr"=agegroup_distribution[ 13 ] )/100
  }
  if ( country != "Australia" ){
    #df_population <- read_excel(  paste0( path_core_functions, "data/eurostat_population.xlsx" ), sheet=3 )
    df_population <- read_fst( paste0( path_core_functions, "data/eurostat_population.fst" ) )
    df_population[ 18, 1 ] <- "Germany" 
    country_row <- which( country == df_population[ ,1 ] )
    country_population <- as.numeric( df_population[ country_row,24 ] )
    agegroup_sizes <- country_population*agegroup_distribution
  }
  return( agegroup_sizes )
}


dem_a <- obtain_demography("Austria")
obtain_demography("Sweden")
c <- contact_matrix(survey = polymod, countries = "Germany", age.limits = c(0, 60), symmetric = TRUE, per.capita = TRUE)
c <- c$matrix.per.capita %*% diag(c(sum(dem_g[1:6]),sum(dem_g[7:10]) )) 
list_surveys()


