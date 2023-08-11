model_SIR_simple = function( params=NULL, data=NULL, country_short_input, agegroup_input, target_input, scenario_tag){
  
    
  data_mock = data %>% 
    filter(country_short == country_short_input, 
           target == target_input, 
           agegroup == agegroup_input)
  warning("Need to fix which data we work with!")
  data_mock %<>% filter(date>"2022-06-01")
  data_mock_fit = data_mock
  # Projection dates
  data_mock_project = data_mock
  data_mock_project$date = data_mock_project$date+365
  
  options(mc.cores = detectCores()-1 )
  set_cmdstan_path(path = NULL)
  mod2 <- cmdstan_model(stan_file='./stan/SIR_simple.stan') # This compiles the script
  
  stan_list = list(
    n_week_fit = nrow(data_mock_fit),
    severe_obs_fit = data_mock_fit$value,
    n_week_project = nrow(data_mock_project),
    pop = 9e6,
    Rnull = params$Rnull,
    rate_infectious = params$rate_infectious
  )
  fit02 <- mod2$sample(
    data = stan_list,
    seed = 12,
    chains = 8,
    parallel_chains = 8,iter_sampling=1500,thin=10,max_treedepth = 15
  )
  
  # Create an output dataframe
  df_out = fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>%
    ungroup() %>%
    left_join(data_mock %>% select(date) %>% mutate(t_vw = 1:n()), by="t_vw") %>%
    mutate(model = "SIR_simple",
           country_short = country_short,
           agegroup = agegroup,
           target = target,
           scenario_tag = scenario_tag,
           prediction_type = "sample"
           ) %>%
    rename(value = .value,
           sample_or_quantile = .draw) %>%
    select(-.chain, -.iteration, -.variable, -t_vw) 
  
    
  if (F){ # plotting to support debugging
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe,
                           pop_infect
    ) %>% 
      mean_qi() -> xp; xp
    
    
    round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
    (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
    (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
    (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    fit02 %>% gather_draws(gen_severe_obs_project[t_vw]) %>% 
      mean_qi() %>% left_join(data_mock_fit %>% mutate(t_vw = 1:n()),by="t_vw") %>%
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      geom_point(aes(y=value),col="black") +
      labs(subtitle = paste("Austria 2022/2023: fit |",
                            "prob_severe:", mprob_severe,"\n",
                            "| S_ini:",mS_ini,
                            "| I_ini:",mI_ini,
                            "| prop inf:",mProp_inf)) -> p_cf0; p_cf0
    
    
  }
  
  return(df_out)
}