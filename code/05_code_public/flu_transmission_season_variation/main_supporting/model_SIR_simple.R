
model_SIR_simple_r0 = function( params=NULL, all_season=NULL, target_input=NULL,pop_country, country_short_input, date_v_fit,season){
  
  # ---- |-Filtering ----
  all_season %>% 
    filter(country_short == country_short_input) %>% 
    select(-typing_sentinel,-typing_nonsentinel,-typing_combined) %>% 
    unnest(inc_iliari) %>% 
    filter(target==params$SIR_simple$target) -> all_season_filtered
  # filter:age groups
  all_season_filtered %>% 
    filter(agegroup==params$SIR_simple$agegroup) -> all_season_filtered
  # filter:time
  all_season_filtered %>% 
    filter( date%in%date_v_fit ) %>% 
    select(country_short,season,date,value) %>% 
    mutate(n=1:n()) -> all_season_fit
  
  if (target_input=="ili_typing_sentinel") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_sentinel)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value)
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing/100)) %>% select(-value_typing)
  }
  if (target_input=="ili_typing_all") {
    xtyping = all_season %>% filter(country_short==country_short_input) %>% 
      unnest(c(typing_combined)) %>% filter(indicator=="positivity") %>% 
      filter( date%in%date_v_fit ) %>% 
      select(date,value_typing=value_add_narm) %>% 
      mutate(value_typing=ifelse(is.nan(value_typing),NA,value_typing ) )
    all_season_fit = all_season_fit %>% left_join(xtyping,by="date") %>% 
      mutate(value=value*(value_typing)) %>% select(-value_typing)
  }
  
  # prepare data for stan fit
  all_season_fit %>% ggplot(aes(date,value)) + geom_line()
  data_mock_fit = all_season_fit
  # Projection dates
  data_mock_project = all_season_fit
  data_mock_project$date = data_mock_project$date+365
  
  #mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = as.integer(replace_na(data_mock_fit$value,0)),
    severe_obs_notna = as.integer(!is.na(data_mock_fit$value)),
    n_week_project = nrow(data_mock_project),
    pop = pop_country,
    Rnull_mu = params$Rnull,
    Rnull =  params$Rnull,
    rate_infectious = params$rate_infectious
  )
  
  if (params$debug==T) {
    est_Rnull = new("precis", .Data = list(c(1.12, 0.08, 
                                 1.00, 1.24, 323.0, 1.0
    )), digits = 2, names = "result", row.names = c("mean", "sd", 
                                                    "5.5%", "94.5%", "n_eff", "Rhat"), .S3Class = "data.frame")
  } else { 
    fit02=rstan::stan(
      file="stan/SIR_simple_nas_freebeta.stan",
      #chains=1 ,thin=1,iter=300, # good for testing model
      chains=14 ,thin=7,iter=560,
      seed=12, cores = getOption("mc.cores", 1L),
      control=list(
        #adapt_delta=0.9,
        #max_treedepth=14
      ),
      data=stan_list
    ) # 46 sec
    est_Rnull_eff = precis(fit02,pars = c("Rnull_eff"))
  }
  # prepared the returned dataframe
  mout = tibble(
    country_short = country_short_input,
    season = season,
    Rnull_eff = est_Rnull_eff$result[1],
    Rnull_eff_Rhat = est_Rnull_eff$result[6]
  )
  
  return(mout)
}