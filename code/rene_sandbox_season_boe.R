# ---- |-Load "infection intensity" data ----
path_data = "https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-truth/ECDC/truth_ECDC-Incident%20Cases.csv"
case_incidenc_EU = read_csv(path_data)
g(case_incidenc_EU)

# filter for minimal working version
case_incidenc_EU = case_incidenc_EU %>% filter(location%in%c("DK",
                                                             "EE","HU","IT",
                                                             "LU","NL","NO","PT","SE"))

# ---- |-fill NA ----
case_incidenc_EU = case_incidenc_EU %>% fill(value,.direction="down")
case_incidenc_EU = case_incidenc_EU %>% group_by(date) %>% summarise(value=sum(value)) %>% 
  mutate(location="EU")

# ---- |-Epiestim to get Rt ----
library(EpiEstim)
n_location = n_distinct(case_incidenc_EU$location)
df_country = list()
for (i in 1:n_location){
  country_i = unique(case_incidenc_EU$location)[i]
  case_incidenc_location = case_incidenc_EU %>% filter_log(location ==country_i) 
  rt <- estimate_R(case_incidenc_location$value, 
                   method="parametric_si",
                   config = make_config(list(
                     mean_si = 1.1, 
                     std_si = 1.5/7) ) )
  filtered_Rt = rt$R$`Mean(R)`
  filtered_Rt[filtered_Rt>2] = NA
  filtered_Rt[filtered_Rt<0.2] = NA
  
  Rt_df = tibble( id=rt$R$t_start,Rt=filtered_Rt )
  
  df_country[[i]] = case_incidenc_location %>% 
    mutate(id=c(1:n())) %>% left_join(Rt_df,by="id")
}

# ---- |-Adding immunity from infection intensity and waning ----


dark_factor = 60 # based on max possible factor (12.1)
NE_factor = 0.8 #
Rt_for_ini_est = 4
R0 = 6
wane_imm_fract_weekly = 0.02
NE = 0.9
NE_and_dark_factor = dark_factor*NE

pop_myexample =120793614

# data needed
aux = list()
aux$pops = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-locations/locations_eu.csv",
                    show_col_types=F)

df_all = bind_rows(df_country) %>% left_join(aux$pops) %>% filter(date>"2023-01-01",date<"2024-02-03")
df_all$population =  pop_myexample
df_all$location_name = "EU"
df_all$Rt_at_zero = 0
df_all$Rt_at_zero[4] = 1
df_all$Rt_at_zero[45] = 2


df_Rt_value_pop = df_all

compute_Rt_without_immunity = function( par, 
                                        df_Rt_value_pop,for_optim=F  ) {
  ## renaming
  NE_and_dark_factor = 80
  wane_imm_fract_weekly_1_log = par[1]
  Imm_2_partial = par[2]
  wane_imm_fract_weekly_2_log = par[3]
  # wane_imm_fract_weekly_3_log = par[4]
  df_all  = df_Rt_value_pop
  
  ## transformations
  wane_imm_fract_weekly_1 = exp(wane_imm_fract_weekly_1_log)
  wane_imm_fract_weekly_2 = exp(wane_imm_fract_weekly_2_log)
  # wane_imm_fract_weekly_3 = exp(wane_imm_fract_weekly_3_log)
  
  location_v = unique(df_all$location_name)
  # country loop
  df_save = list()
  for (i in location_v) { # i = location_v[1]
    df_i = df_all %>% filter(location_name==i) 
    # normalise incidence
    df_i$value = df_i$value * (0.01/sum( (df_i$value/df_i$population) , na.rm=T  )  )
    # if(sum(df_i$value/df_i$population) != 0.01) print(paste("Check normalising inc for",i))
    # initial state
    Rt_ini = mean(df_i$Rt[1:Rt_for_ini_est])
    Imm_ini = (R0 - Rt_ini)/R0
    # 
    df_i$Imm_1 = NA; df_i$Imm_1[1] = Imm_ini # full immunity
    df_i$Imm_2 = NA; df_i$Imm_2[1] = 0 # partial imminty
    df_i$Imm = NA; df_i$Imm[1] = df_i$Imm_1[1] + df_i$Imm_2[1]*Imm_2_partial # all immunity combined
    #
    df_i$waning_1 = NA; df_i$waning_1[1] = 0
    df_i$waning_2 = NA; df_i$waning_2[1] = 0
    df_i$new_Imm = NA; df_i$new_Imm[1] = 0
    df_i$Rt_without_immunity = NA; 
    df_i$Rt_without_immunity[1] = df_i$Rt[1] * 1/(1-Imm_ini)
    #
    # time loop (through the weeks of a given location)
    for (t in 2:nrow(df_i) ) { # t = 2
      # immunity waning (allowing different immunity compartments)
      # Imm_1 is the full immunity compartment
      # Imm_2 3, ... are partial immunity compartments
      waning_1 = df_i$Imm_1[t-1]*wane_imm_fract_weekly_1 
      waning_2 = df_i$Imm_2[t-1]*wane_imm_fract_weekly_2 
      df_i$waning_1[t] = waning_1
      df_i$waning_2[t] = waning_2
      # accumulatoin of new immunity (proportional to: value = "infection intensity",dark_factor,NE)
      new_Imm = df_i$value[t-1]/df_i$population[t-1] * NE_and_dark_factor
      df_i$new_Imm[t] = new_Imm
      # apply changes to Immunity compartments
      df_i$Imm_1[t] = df_i$Imm_1[t-1] - waning_1 + new_Imm
      df_i$Imm_2[t] = df_i$Imm_2[t-1] - waning_2 + waning_1
      if (df_i$Imm_1[t]<0) print("warning1")
      if (df_i$Imm_2[t]<0) print("warning2")
      # recompute combined "effective" immunity
      df_i$Imm[t] = df_i$Imm_1[t] + df_i$Imm_2[t]*Imm_2_partial
      # compute new Rt without immunity (Rt AS IF immunity would have vanished)
      df_i$Rt_without_immunity[t] = df_i$Rt[t] * 1/(1-df_i$Imm[t])
      
      ## scenarios 
      # df_i$s_1[t] =  df_i$s_1[t-1] + new_Imm
    }
    # target allows multiple contraints
    # regularise: immunity last year = immunity 1 year later
    Imm_diff = df_i$Imm[df_i$Rt_at_zero==1] - df_i$Imm[df_i$Rt_at_zero==2]
    target_Imm_diff = -dnorm(x=Imm_diff,mean = 0,sd=0.05, log = T)
    target_Imm_diff_min = -dnorm(x=0,mean = 0,sd=0.05, log = T)
    # regularise: R0 based on Rt and R0
    max_Rt_without_immunity = max(df_i$Rt_without_immunity,na.rm=T)
    target_Rt_without_immunity = -dnorm(x=max_Rt_without_immunity,mean = R0,sd=1, log = T)
    target_Rt_without_immunity_min = -dnorm(x=R0,mean = R0,sd=1, log = T)
    # regulatise; max immunity in the population
    
    #
    target = target_Imm_diff + target_Rt_without_immunity
    min_target = target_Imm_diff_min + target_Rt_without_immunity_min # -1.157855
    # print(paste("min target:",min_target %>% round(2)))
    # print(paste("Max new_Imm",max(df_i$new_Imm) %>% round(3),i))
    # print(paste("Average waning",mean(df_i$waning) %>% round(3),i))
    df_save[[i]] = df_i
  }
  df_return = bind_rows(df_save)
  if (for_optim==T) return(target)
  if (for_optim==F) return(df_return)
}


# ---- |-Find waning function and NE_factor that work ----
# initial para guess
# pars=c(log(0.02),0.5 , log(0.02))
# ( opt_out = optim(par=pars,method="L-BFGS-B",
#                   lower=c(log(0.005),0.1,log(0.005)),upper=c(log(0.1),0.99,log(0.1) ) , 
#                   fn=compute_Rt_without_immunity,df_Rt_value_pop=df_Rt_value_pop,for_optim=T) )
# opt_out$par
# df_Rt_without_imm = compute_Rt_without_immunity( par=opt_out$par , 
#                                                  df_Rt_value_pop , 
#                                                  for_optim=F)

df_min = df_grid_all %>% filter(target==min(target)) 
df_Rt_without_imm = compute_Rt_without_immunity(par=c(df_min$wane_imm_fract_weekly_1_log,
                                                      df_min$Imm_2_partial,
                                                      df_min$wane_imm_fract_weekly_2_log) , 
                                                df_Rt_value_pop , 
                                                for_optim=F)



# df_all %<>% group_by(location) %>% 
#   filter(date>"2023-01-01",date<"2023-11-15") %>% 
#   mutate(
#     # normalise incidence
#     value = value *(0.01/sum( (value/population) , na.rm=T) ),
#     # get Rt correction based on accumumulating immunity
#     value_cumsum = cumsum(value),
#     fract_case_inverse = 1/(value_cumsum/population),
#     fract_infect = (value_cumsum*dark_factor)/population,
#     fract_immunised = (value_cumsum*dark_factor*NE_factor)/population,
#     Rt_reduction = (1-fract_immunised),
#     Rt_reduction_inverse =1/(Rt_reduction),
#     Rt_without_immunity_old = Rt*Rt_reduction_inverse,
#     # approach 2 (using value instead of cumsum(value) but then cumprod)
#     Rt_reduction_inverse_2 = cumprod(1/(1-(value*dark_factor*NE_factor)/population)),
#     Rt_without_immunity = Rt*Rt_reduction_inverse_2,
#   ) #%>% summarise(check=max(Rt_reduction_inverse)) %>% print(n=30)
# df_all %>% summarise(check=max(Rt_without_immunity)) %>% print(n=Inf)
# # ---- |-Look at data ----


# ---- |-Prepare df for stan ----
df_stan = df_Rt_without_imm %>% ungroup() %>% 
  filter(date>"2023-01-01",date<"2023-11-15") %>% filter(!is.na(Rt)) %>% 
  mutate(season_t=as.numeric(date-ymd("2022-12-31")) ) %>% 
  mutate(id=c(1:n()))

# ---- |-Create stan list ----
stan_list = list(
  # model input related to data-points, country structure, etc
  n = nrow(df_stan),
  n_location = n_distinct(df_stan$location),
  location_id = fct_inorder(df_stan$location) %>% as.numeric(),
  location_id_raw = fct_inorder(df_stan$location) %>% levels() %>% enframe,
  # model input related to seasonality and outcome
  season_t = df_stan$season_t,
  y = df_stan$Rt_without_immunity,
  # model input related to prior information
  prior = 2
)

# ---- |-Create stan list ----
fit_season=rstan::stan(
  file="./stan/season00_start.stan",
  chains=8 ,thin=8,iter=400,# 1200,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    # max_treedepth=14
  ),
  data=stan_list
) # 
save(fit_season,stan_list,file="../Big data/fit_season.Rdata")
load(file="../Big data/fit_season.Rdata")
fit = fit_season

library(rethinking)
precis(fit,pars=c("season_fold","season_amp"),depth=2) # 40-50% more beta in winter compared to summer
precis(fit,pars=c("season_fold_log_mu","season_fold_log_sd","y_base","sigma"),depth=2)

# ---- |-Predictive checks ----
df_pp = fit %>% gather_draws(mu[id]) %>% filter(.draw%in%c(1:50)) %>% 
  left_join(df_stan , by="id")

df_stan %>% 
  ggplot(aes(x=date,y=Rt_without_immunity)) + 
  geom_line() + geom_hline(yintercept = 1) + 
  geom_line(data=df_pp,aes(x=date,y=.value,group=.draw),alpha=0.2,col="blue" ) +
  facet_wrap(~location,scales="free_y") + 
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())

# 1.47 fold increase from summer to winter
# =47% increase from summer to winter
# =32% decrease from winter to summer
# 42.1% reduction from winter to summer (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010435)
