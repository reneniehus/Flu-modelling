# ---- |-Setup: cmdstanr ----
library(cmdstanr)
options(mc.cores = detectCores()-1 )
set_cmdstan_path(path = NULL)
mod <- cmdstan_model(stan_file='./stan/immu_dyn_00.stan')

# ---- |-Simulate: "Australian wave" ----
wave_settings_SIR_Rnull = list(
  I_ini=0.0001,
  S_ini= 0.8,
  Rnull = 2.5,
  prop_severe = 0.02,
  sigma_m = log(50)*0.05
)
wave_df_w = sim_wave_SIR_Rnull_w(wave_settings_SIR_Rnull,"test")
wave_df_w = sim_wave_SIR_Rnull_w(wave_settings_SIR_Rnull,"test")

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
fit01 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 8,
  parallel_chains = 8,iter_sampling=5000,thin=5,max_treedepth = 15
)
# get_variables(fit01)
precis(fit01,pars = c("I_ini","prop_severe","S_ini"))