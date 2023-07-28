# ---- |-Setup: cmdstanr ----
library(cmdstanr)
options(mc.cores = detectCores()-1 )
set_cmdstan_path(path = NULL)
#mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart.stan') # This compiles the script
mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v2.stan') # This compiles the script


# Load EU flu data from past seasons
country <- "BE"
source <- "ILI"
connectionStringSource <- "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDM;Trusted_Connection=TRUE;"
channelSource = odbcDriverConnect(connectionStringSource)
sqlCode <-
  paste(
    "EXEC [dbo].[spMEMGetInflClin] @Region = N'",
    country,
    "', @Numerator = N'",
    source,
    "'",
    sep = ''
  )

df <- sqlQuery(channelSource, sqlCode) %>% as.data.frame()
df$t_vw <- 1:nrow(df)
inputData <- df %>% 
  pivot_longer(cols = colnames(df)[2:(length(colnames(df))-1)], names_to = "Season", values_to = "value")

inputData %>% 
  ggplot(aes(x=t_vw, y=value, color=Season)) + geom_line() + scale_y_log10()

inputData$value[is.na(inputData$value)] <- 0
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Australia + EU starting off ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---- |-Plot: Australian wave + start of EU wave  ----
df_to_fit <- inputData %>% 
  filter(Season=="S2014-2015") %>% 
  filter(value>100)
df_to_fit$t_vw = c(1:nrow(df_to_fit))
df_to_fit %>%  ggplot(aes(x=t_vw,y= value )) + 
  geom_line() #+ scale_y_log10()
df_to_fit$value_infl = df_to_fit$value


# wave_start = df_aus %>% mutate(value_infl=value_infl,
#                   value_infl_cumsum = cumsum(value_infl) ) %>% 
#   filter(value_infl_cumsum<2000)
wave_start0 = inputData %>% 
  filter(Season=="S2017-2018") %>%
  mutate(value_infl=value,
         value_infl_cumsum = cumsum(value_infl) ) %>%
  filter(value_infl > 50)

wave_start = wave_start0 %>% 
  filter(value_infl_cumsum<1000)
wave_start$t_vw = c(1:nrow(wave_start))
wave_start0$t_vw = c(1:nrow(wave_start0))
wave_start %>%  ggplot(aes(x=t_vw,y= value )) + 
  geom_point() #+ scale_y_log10()


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
### AFit: Australian wave + start of EU wave  ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pop_fit = 11.6e6 # pop_australia
pop_target = 11.6e6 #pop_australia
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
# get_variables(fit02)
precis(fit02,pars = c("SIR_ini","I_ini_start","prop_severe","prop_severe_start",
                      "pop_infect","pop_infect_start"), depth=3) # includes rhat diagnostics
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

round(xp[xp$.variable=="prop_severe_start"&xp$scenario==2,".value"],4) -> mprob_severe
(xp[xp$.variable=="I_ini_start"&xp$scenario==2,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==2) %>% 
  mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_line(aes(y=.value)) +
  #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
  geom_point(data=wave_start0,aes(x=t_vw, y=value_infl )) +
  geom_point(data=wave_start,aes(y=value_infl),col="darkred") + 
  #scale_y_log10() + 
  labs(subtitle = paste("EU scenario: more transmissible |","\n",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini) ) -> p_cf2; p_cf2

round(xp[xp$.variable=="prop_severe_start"&xp$scenario==3,".value"],4) -> mprob_severe
(xp[xp$.variable=="I_ini_start"&xp$scenario==3,".value"]) %>% logit() %>% round(1) -> mI_ini
fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==3) %>% 
  mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
  ggplot(aes(x=t_vw)) + 
  geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
  geom_line(aes(y=.value)) +
  #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
  geom_point(data=wave_start0,aes(x=t_vw, y=value_infl )) +
  geom_point(data=wave_start,aes(y=value_infl),col="darkred")+
  #scale_y_log10() + 
  labs(subtitle = paste("EU scenario: more vaccine |","\n",
                        "prob_severe:", mprob_severe,
                        "| I_ini:",mI_ini) ) -> p_cf3; p_cf3

mexplain = c("Fitting Australian wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
             prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
More transmissible: fit with higher beta
More vaccination: fit with lower S_ini,
grey bands showing 80% credible intervals")
#( (p_cf0 / (p_cf1 + p_cf2 + p_cf3) ) + plot_annotation(caption = mexplain) )-> ppatch;ppatch
( ( (p_cf0+p_cf1) / ( p_cf2 + p_cf3) ) + plot_annotation(caption = mexplain) )-> ppatch;ppatch
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






# ---- |-Fit: real data  ---- # FIXME: can delete
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
