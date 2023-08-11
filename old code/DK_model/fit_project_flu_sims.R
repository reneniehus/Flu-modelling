source("setup01.R")
source("DK_model/bayesian_functions.R")

# ---- |-Setup: cmdstanr (we are using cmdstanr instead of rstan because it has more features) ----
library(cmdstanr)
options(mc.cores = detectCores()-1 ) # detect number of cores for parallel chains
set_cmdstan_path(path = NULL)
mod <- cmdstan_model(stan_file='stan/immu_dyn_00.stan') # compiles the stan model into an .exe

# ---- |-Simulate: "Australian wave" ----
wave_settings_SIR_Rnull = list(
  I_ini=0.0001,
  S_ini= 0.8,
  Rnull = 2.5,
  prop_severe = 0.02,
  sigma_m = log(50)*0.05
)
wave_df = sim_wave_SIR_Rnull(wave_settings_SIR_Rnull,"test")
wave_df_w = sim_wave_SIR_Rnull_w(wave_settings_SIR_Rnull,"test")
wave_df %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean)) + scale_y_log10() + 
  labs(caption="test")
wave_df_w %>% ggplot(aes(x=i_day,y=severe_obs_weekly)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean_weekly)) + #scale_y_log10() + 
  labs(caption="test")

# ---- |-Fit: "Australian wave"  ----
df_hosp_sel = wave_df_w
stan_list = list(
  n_week = nrow(df_hosp_sel),
  severe_obs = (df_hosp_sel$severe_obs_weekly) %>% replace_na(replace = 0),
  pop=df_hosp_sel$pop[1]
)
# here the model is being fit:
# seed is the seed for pseudorandom numbers and allow to reproduce exactly the same posterior
# chains gives number of independently run chains (we want >1, and may as well use our processor kernels)
# iter_sampling gives the samples that we get after warm up phase (which is set as a default if not specified)
# thin=5 means we only keep every 5th sample to have less stuff in memory
# max_treedepth is a parameter of the HMC search and requires adaptation if there are many divergent transisions
fit01 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 8,
  parallel_chains = 8,iter_sampling=5000,thin=5,max_treedepth = 15
)
# get_variables(fit01)
precis(fit01,pars = c("I_ini","prop_severe","S_ini"))
fit01 %>% gather_draws(I_ini_logit,prop_severe_logit,S_ini_logit,
                             I_ini_logit_prior,prop_severe_logit_prior,S_ini_logit_prior) %>% 
  mean_qi() %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = c(logit(wave_settings_SIR_Rnull$I_ini),
                            logit(wave_settings_SIR_Rnull$prop_severe),
                            logit(wave_settings_SIR_Rnull$S_ini)), linetype = 'dotted')
# generated dynamics
fit01 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==1) %>% 
  mean_qi() %>% left_join(wave_df_w,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_point(aes(y=severe_obs_weekly))+labs(subtitle = "Fake Australia: fitted") -> p_cf0

fit01 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==2) %>% 
  mean_qi() %>% left_join(wave_df_w,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_point(aes(y=severe_obs_weekly)) +labs(subtitle = "More transmissible") -> p_cf1

fit01 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==3) %>% 
  mean_qi() %>% left_join(wave_df_w,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_point(aes(y=severe_obs_weekly)) +labs(subtitle = "More vaccination") -> p_cf2
p_cf0 + p_cf1 + p_cf2


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Australia + EU starting off ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---- |-Simulate: Australian wave + start of EU wave  ----
mod2 <- cmdstan_model(stan_file='stan/immu_dyn_01_eustart.stan')
#
wave_settings_SIR_Rnull = list(
  I_ini=0.000001,
  S_ini= 0.7,
  Rnull = 2.0,
  prop_severe = 0.02,
  sigma_m = log(50)*0.05
)
wave_df_w = sim_wave_SIR_Rnull_w(wave_settings_SIR_Rnull,"EU")
wave_full = wave_df_w
wave_start = wave_df_w %>% 
  mutate(severe_obs_weekly_cumsum=cumsum(severe_obs_weekly)) %>% 
  filter(severe_obs_weekly_cumsum<200)
wave_full %>% ggplot(aes(x=i_day,y=severe_obs_weekly)) + 
  geom_point(col=NA)+
  geom_point(data=wave_start,col="darkred")+
  labs(subtitle ="So-far-observed EU wave (red)") -> p_eu
# Australia
wave_settings_SIR_Rnull$I_ini=0.001
wave_australia = sim_wave_SIR_Rnull_w(wave_settings_SIR_Rnull,"Australia")
wave_australia %>% ggplot(aes(x=i_day,y=severe_obs_weekly)) + geom_point(alpha=0.2) + 
  geom_point(alpha=1)+
  labs(subtitle="Australia wave") ->p_aust
mexplain = c("Simulated data for 2022 Australia wave and a beginning autumn 2022 EU wave:
             assuming single, well-mixed population, 
fixed values for pandemic flu R0, and serial interval 
Time-discrete (daily) simulation with Australia settings: 
I_ini=0.1%, S_ini=0.7, prop_severe=0.02 and little noise,
EU settings: I_ini=0.0001%, otherwise identical parameters")
( p_aust / p_eu ) + plot_annotation(caption = mexplain)

# ---- |-Fit: Australian wave + start of EU wave  ----
stan_list = list(
  n_week_full = nrow(wave_australia),
  severe_obs_full = (wave_australia$severe_obs_weekly) %>% replace_na(replace = 0),
  n_week_start = nrow(wave_start),
  severe_obs_start = (wave_start$severe_obs_weekly) %>% replace_na(replace = 0),
  pop_full=wave_full$pop[1],
  pop_start=wave_start$pop[1]
)
fit02 <- mod2$sample(
  data = stan_list,
  seed = 12,
  chains = 8,
  parallel_chains = 8,iter_sampling=10000,thin=10,max_treedepth = 15
)

# ---- |-Plot: Parameters  ----
# get_variables(fit02)
precis(fit02,pars = c("SIR_ini","I_ini_start","prop_severe","prop_severe_start")) # includes rhat diagnostics

fit02 %>% gather_draws(SIR_ini[state],I_ini_start[scenario],prop_severe,prop_severe_start[scenario]) %>% 
  mean_qi() -> xp; xp


fit02 %>% gather_draws(I_ini_logit,prop_severe_logit,S_ini_logit,
                       I_ini_logit_prior,prop_severe_logit_prior,S_ini_logit_prior) %>% 
  mean_qi() %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = c(logit(wave_settings_SIR_Rnull$I_ini),
                            logit(wave_settings_SIR_Rnull$prop_severe),
                            logit(wave_settings_SIR_Rnull$S_ini)), linetype = 'dotted')

# ---- |-Plot: Scenario projections  ----
round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
(xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs_full[t_vw]) %>% 
  mean_qi() %>% left_join(wave_australia,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_point(aes(y=severe_obs_weekly),col="black")+
  labs(subtitle = paste("Australia: fit |",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini)) -> p_cf0; p_cf0

round(xp[xp$.variable=="prop_severe_start"&xp$scenario==1,".value"],3) -> mprob_severe
(xp[xp$.variable=="I_ini_start"&xp$scenario==1,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==1) %>% 
  mean_qi() %>% left_join(wave_full,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
  geom_point(data=wave_start,aes(y=severe_obs_weekly),col="darkred")+
  labs(subtitle = paste("EU scenario: baseline |",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini) ) -> p_cf1; p_cf1

round(xp[xp$.variable=="prop_severe_start"&xp$scenario==2,".value"],4) -> mprob_severe
(xp[xp$.variable=="I_ini_start"&xp$scenario==2,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==2) %>% 
  mean_qi() %>% left_join(wave_full,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
  geom_point(data=wave_start,aes(y=severe_obs_weekly),col="darkred")+
  labs(subtitle = paste("EU scenario: more transmissible |",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini) ) -> p_cf2; p_cf2

round(xp[xp$.variable=="prop_severe_start"&xp$scenario==3,".value"],4) -> mprob_severe
(xp[xp$.variable=="I_ini_start"&xp$scenario==3,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==3) %>% 
  mean_qi() %>% left_join(wave_full,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
  geom_point(data=wave_start,aes(y=severe_obs_weekly),col="darkred")+
  labs(subtitle = paste("EU scenario: more vaccine |",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini) ) -> p_cf3; p_cf3

mexplain = c("Fitting Australian wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
             prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
More transmissible: fit with higher beta
More vaccination: fit with lower S_ini,
grey bands showing 80% credible intervals")
(p_cf0 / (p_cf1 + p_cf2 + p_cf3) ) + plot_annotation(caption = mexplain)



# ---- |-Fit: real data  ----
stan_list = list(
  n_week = nrow(df_stan),
  severe_obs = as.numeric(df_stan$INF_ALL) %>% replace_na(replace = 0),
  pop=pop_denmark
)
fit_den <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 8,
  parallel_chains = 8,iter_sampling=5000,thin=5,max_treedepth = 15
)
