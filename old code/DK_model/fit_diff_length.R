df_hosp = sim_wave_SIR(NULL)

# ---- |-Plotting: different data subsets ----
df_hosp$t_v %>% quantile(prob=c(0.25,0.5,0.75))
df_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.1)+
  geom_line(aes(y=severe_mean),col="blue") +
  geom_point(data=. %>% filter(t_v<61))
df_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.1)+
  geom_line(aes(y=severe_mean),col="blue") +
  geom_point(data=. %>% filter(t_v<31))
df_hosp %>% ggplot(aes(x=t_v,y=severe_obs)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.1)+
  geom_line(aes(y=severe_mean),col="blue") +
  geom_point(data=. %>% filter(t_v<91))

# myinitofit <- lapply(1:4, function(id) list( beta_logit=-0.847,
#                                              prop_hosp_logit=-3.891,
#                                              S_ini_logit=1.38,
#                                              simga=0.02) )

# ---- |-Setuo: cmdstanr ----
library(cmdstanr)
options(mc.cores = detectCores()-1 )
set_cmdstan_path(path = NULL)
mod <- cmdstan_model(stan_file='./stan/immu_dyn_00.stan')

# ---- |-Fits: source of the bias ----
df_hosp_sel = df_hosp
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs = (df_hosp_sel$severe_obs) %>% replace_na(replace = 0),
  pop=df_hosp_sel$pop[1]
)
fit100_perf <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 8,
  parallel_chains = 8,iter_sampling=5000,thin=5,max_treedepth = 15
)
# get_variables(fit100_perf)

precis(fit100_perf, pars = c("beta","prop_severe","S_ini")
       )
#save(fit100_perf,file = "./stanfit/immu_dyn.Rdata")
fit100_perf %>% gather_draws(beta_logit,prop_severe_logit,S_ini_logit,
                             beta_logit_prior,prop_severe_logit_prior,S_ini_logit_prior) %>% 
  mean_qi() %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = c(logit(beta),logit(prop_hosp),logit(S_ini)), linetype = 'dotted')

fit100_perf %>% gather_draws(severe_obs_mean[time]) %>% 
  mean_qi() -> df_gen_hosp

tibble(x=c(1:nrow(df_hosp_sel)),
       y_data=df_hosp_sel$severe_obs,
       y_model=df_gen_hosp$.value) %>% 
  ggplot(aes(x=x)) +
  geom_line(aes(y=(y_data))) + geom_point(aes(y=(y_model)))
plot(df_gen_hosp$.value, (df_hosp_sel$delta_severe*pop)  )

# look at how model's SEIR looks like
fit100_perf %>% spread_draws(S[time],I[time],E[time],R[time]) %>% 
  filter(.chain==max(.chain),.draw==max(.draw))
df_hosp_sel

# look at correlation of draws
fit100_perf %>% spread_draws(beta_logit,prop_severe_logit,S_ini_logit) %>% 
  select(-c(1:3)) %>% pairs()


df_gen_hosp$.value %>% plot()
(df_hosp_sel$delta_severe*pop) %>% plot()
# ---- |-Fits: to different subsets of the data ----
# 100% no noise
df_hosp_sel = df_hosp
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs_log = log(df_hosp_sel$delta_severe*pop),
  pop=df_hosp_sel$pop[1]
)
fit100_perf <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 4,
  parallel_chains = 4
)

# 100%
df_hosp_sel = df_hosp
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs_log = log(df_hosp_sel$severe_obs),
  pop=df_hosp_sel$pop[1]
)
fit100 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 4,
  parallel_chains = 4
)
# 75%
df_hosp_sel = df_hosp %>% filter(t_v<91)
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs_log = log(df_hosp_sel$severe_obs),
  pop=df_hosp_sel$pop[1]
)
fit075 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 4,
  parallel_chains = 4
)
# 50%
df_hosp_sel = df_hosp %>% filter(t_v<61)
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs_log = log(df_hosp_sel$severe_obs),
  pop=df_hosp_sel$pop[1]
)
fit050 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 4,
  parallel_chains = 4
)
# 25%
df_hosp_sel = df_hosp %>% filter(t_v<31)
stan_list = list(
  n = nrow(df_hosp_sel),
  severe_obs_log = log(df_hosp_sel$severe_obs),
  pop=df_hosp_sel$pop[1]
)
fit025 <- mod$sample(
  data = stan_list,
  seed = 12,
  chains = 4,
  parallel_chains = 4
)

# beta = 0.3; logit(beta)
# prop_hosp = 0.02; logit(prop_hosp)
# S_ini = 0.80; logit(S_ini)
# sigma_m = 0.02

precis(fit100_perf)
pairs(fit100)

plot_this = function(df) {
  p = df %>% 
    gather_draws(beta_logit,prop_hosp_logit,S_ini_logit) %>% 
    mean_qi() %>% 
    ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
    geom_pointrange() +
    geom_vline(xintercept = c(logit(beta),logit(prop_hosp),logit(S_ini)), linetype = 'dotted')
  return(p)
}


fit100_perf %>% gather_draws(beta_logit,prop_hosp_logit,S_ini_logit,
                             beta_logit_prior,prop_hosp_logit_prior,S_ini_logit_prior) %>% 
  mean_qi() %>% 
  ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
  geom_pointrange() +
  geom_vline(xintercept = c(logit(beta),logit(prop_hosp),logit(S_ini)), linetype = 'dotted')


fit100%>% plot_this()
fit075%>% plot_this()
fit050%>% plot_this()
fit025%>% plot_this()






