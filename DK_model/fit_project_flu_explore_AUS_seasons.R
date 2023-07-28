# ---- |-Setup: cmdstanr ----
library(cmdstanr)
options(mc.cores = detectCores()-1 )
set_cmdstan_path(path = NULL)
#mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart.stan') # This compiles the script
mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v4.stan') # This compiles the script


for (y in 0:9){
  
  yr_start = paste0("201",y,"-03-01")
  yr_end = paste0("201",y,"-12-01")
  aus_ili = read_csv(file="./data/data_AUS_ILI.csv")
  aus_ili %>% filter(date> yr_start, date<yr_end) %>% ggplot(aes(x=date,y=value)) + geom_line()
  aus_ili %>% filter(date> yr_start, date<yr_end) -> df_aus
  
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Australia + EU starting off ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # ---- |-Plot: Australian wave + start of EU wave  ----
  df_to_fit <- df_aus
  df_to_fit$t_vw = c(1:nrow(df_to_fit))
  df_to_fit %>%  ggplot(aes(x=t_vw,y= value )) + 
    geom_line() #+ scale_y_log10()
  df_to_fit$value_infl = df_to_fit$value
  
  
  # wave_start = df_aus %>% mutate(value_infl=value_infl,
  #                   value_infl_cumsum = cumsum(value_infl) ) %>% 
  #   filter(value_infl_cumsum<2000)
  wave_start = df_to_fit
  
  
  
  
  df_to_fit %>% ggplot(aes(x=t_vw,y=value_infl)) +
    geom_point(col=NA)+
    geom_point(data=wave_start,col="darkred")+
    labs(subtitle ="So-far-observed EU wave (red)") -> p_eu;p_eu
  
  # wave_start %>% ggplot(aes(x=1:length(value_infl),y=value_infl)) + 
  #   geom_point()+
  #   #geom_point(data=wave_start,col="darkred")+
  #   labs(subtitle ="So-far-observed EU wave (red)") -> p_eu;p_eu
  
  
  # Australia
  df_to_fit %>% ggplot(aes(x=t_vw,y=value_infl)) + geom_point(alpha=0.2) + 
    geom_point(alpha=1)+
    labs(subtitle="Australia wave") ->p_aust ; p_aust
  mexplain = c("Real ILI data from Australia, multiplied with 100 for poisson model")
  ( p_aust / p_eu ) + plot_annotation(caption = mexplain)
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Fit: Australian wave + start of EU wave  ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pop_fit = 26e6 # pop_australia
  pop_target = 26e6 #pop_australia
  df_for_fit = df_to_fit #df_aus
  stan_list = list(
    n_week_full = nrow(df_for_fit),
    severe_obs_full = (df_for_fit$value_infl) %>% replace_na(replace = 0) %>% round(0),
    n_week_start = nrow(wave_start),
    severe_obs_start = (wave_start$value_infl) %>% replace_na(replace = 0) %>% round(0),
    pop_full=pop_fit,
    pop_start=pop_target
  )
  fit02 <- mod2$sample(
    data = stan_list,
    seed = 12,
    chains = 8,
    parallel_chains = 8,iter_sampling=1500,thin=10,max_treedepth = 15
  )
  
  # ---- |-Plot: Parameters  ----
  # create table of parameters
  fit02 %>% gather_draws(SIR_ini[state],
                         I_ini_start[scenario],
                         prop_severe,
                         pop_infect,
                         #prop_severe_start,
                         pop_infect_start[scenario]) %>% 
    mean_qi() -> xp; xp
  
  
  fit02 %>% gather_draws(I_ini_logit,
                         prop_severe_logit,
                         S_ini_logit,
                         I_ini_logit_prior,
                         prop_severe_logit_prior,
                         S_ini_logit_prior) %>% 
    mean_qi() %>% 
    ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
    geom_pointrange() +
    geom_vline(xintercept = c(logit(wave_settings_SIR_Rnull$I_ini),
                              logit(wave_settings_SIR_Rnull$prop_severe),
                              logit(wave_settings_SIR_Rnull$S_ini)), linetype = 'dotted')
  
  # ---- |-Plot: Scenario projections  ----
  round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
  (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
  (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
  (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
  fit02 %>% gather_draws(gen_severe_obs_full[t_vw]) %>% 
    mean_qi() %>% left_join(df_for_fit,by="t_vw") %>% 
    ggplot(aes(x=t_vw)) + 
    geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
    geom_line(aes(y=.value)) +
    geom_point(aes(y=value_infl),col="black")+
    labs(subtitle = paste("Australia: fit |",
                          "prob_severe:", mprob_severe,"\n",
                          "| S_ini:",mS_ini,
                          "| I_ini:",mI_ini,
                          "| prop inf:",mProp_inf)) -> p_cf0; p_cf0
  
  round(xp[xp$.variable=="prop_severe_start"&xp$scenario==1,".value"],3) -> mprob_severe
  (xp[xp$.variable=="I_ini_start"&xp$scenario==1,".value"]) %>% logit() %>% round(1) -> mI_ini
  fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==1) %>% 
    mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
    ggplot(aes(x=t_vw)) + 
    geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
    geom_line(aes(y=.value)) +
    #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
    geom_point(data=wave_start0,aes(x=t_vw, y=value_infl )) +
    geom_point(data=wave_start,aes(y=value_infl ),col="darkred") +
    #scale_y_log10() + 
    labs(subtitle = paste("EU scenario: baseline |","\n",
                          "prob_severe:", mprob_severe,
                          "| I_ini:",mI_ini) ) -> p_cf1; p_cf1
  
  
  
  mexplain = c("Fitting Australian wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
             prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
More transmissible: fit with higher beta
More vaccination: fit with lower S_ini,
grey bands showing 80% credible intervals")
  #( (p_cf0 / (p_cf1 + p_cf2 + p_cf3) ) + plot_annotation(caption = mexplain) )-> ppatch;ppatch
  ( ( (p_cf0+p_cf1) ) + plot_annotation(caption = mexplain) )-> ppatch;ppatch
  #ggsave_as(ppatch,"real_aus_eu",30,30)
  
  here()
  
  
  # Compare cumulative between data and model
  # Model
  fit02 %>% 
    gather_draws(gen_severe_obs[scen,t_vw]) %>% 
    filter(scen==1) %>% 
    mean_qi() %>%
    pull(.value) %>% 
    sum()
  
  # Data
  wave_start0$value %>% sum()
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Save parameters from the historical fit ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  params_hist <- xp
  S_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==1,".value"]
  I_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==2,".value"]
  R_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==3,".value"]
  prop_severe_hist <- params_hist[params_hist$.variable=="prop_severe",".value"]
  
  
  
  df_tmp <- as_tibble(bind_cols(S_init_hist, I_init_hist, R_init_hist, prop_severe_hist, Y)) %>% 
    rename(S_init = .value...1, I_init = .value...2, R_init = .value...3, prop_severe = .value...4, season = ...5)
  
  df_param <- bind_rows(df_param, df_tmp)
  
  
  
}
