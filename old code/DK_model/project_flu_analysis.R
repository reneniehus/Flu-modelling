# Script to do list:
# - select one country to do the analysis on
# - check MEM epidemic threshold - use only data above it (or perhaps few data point before)
# - fit one historic wave; obtain prop_severe, S_init, I_init, prop_inf
# - compare S_init to AUS 2022 wave - decide which one to use
# - check how well it fits (fitting only I_init?) to another pre-COVID-19 historic wave; check cumulative cases 
# - check different scenarios; get relative decrease in burden for different scenarios
# - sensitivity analysis; vary S_init, prop_severe

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Load libraries and data ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ---- |-Load libraries  ----
source("./code/setup01.R")
source("./code/Flu/DK_model/bayesian_functions.R")
library(cmdstanr)
#install.packages("RODBC")
library(RODBC)
path_core_functions <- "../COVID-19-input/"
source( paste0( path_core_functions, "core_functions/function_obtain_demography.R" ))


# ---- |-Load AUS data  ----
aus_ili = read_csv(file="./data/data_AUS_ILI.csv")
aus_ili %>% 
  #filter(date>"2022-03-01") %>% 
  ggplot(aes(x=date,y=value)) + geom_line()
aus_ili %>% filter(date>"2022-03-01") -> df_aus

# ---- |-Load EU data  ----
#country_short <- "EL"
#country_long <- "Greece"

# COUNTRIES_LONG <- c("Slovenia", "Austria","Latvia","France", "Greece")
# COUNTRIES_SHORT <- c("SI", "AT", "LV", "FR", "EL")
COUNTRIES_LONG <- c("Greece", "Netherlands")
COUNTRIES_SHORT <- c("EL", "NL")
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Spectral")


for (j in 1:2){
  country_short <- COUNTRIES_SHORT[j]
  country_long <- COUNTRIES_LONG[j]
  
  #
  # Old way
  if (F){
    source <- "ILI"
    connectionStringSource <- "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDM;Trusted_Connection=TRUE;"
    channelSource = odbcDriverConnect(connectionStringSource)
    sqlCode <-
      paste(
        "EXEC [dbo].[spMEMGetInflClin] @Region = N'",
        country_short,
        "', @Numerator = N'",
        source,
        "'",
        sep = ''
      )
    df <- sqlQuery(channelSource, sqlCode) %>% as.data.frame()
    df$t_vw <- 1:nrow(df)
    inputData <- df %>% 
      pivot_longer(cols = colnames(df)[2:(length(colnames(df))-1)], names_to = "Season", values_to = "value")
    
  }
  
  # New way
  INFL_Snapshot <- "Driver=SQL Server;Server=nsql3;Database=INFL"
  INFL.connection <- odbcDriverConnect(INFL_Snapshot)
  sql.code <- "SELECT * FROM [INFL].[clean].[INFLCLIN_Haggregated]"
  df0 <- sqlQuery(INFL.connection, sql.code)
  df <- df0 %>% filter(ReportingCountry == country_short) %>%
    select(CountryName, ReportingCountry, Season, TimeCode, DateUsedForStatisticsISO, DateUsedForStatisticsWeek, isSyndromic,
           ILI_DenominatorNumberOfCases, ILINumberOfCases, PercPositive) %>%
    group_by(Season) %>%
    mutate(t_vw = 1:length(ILINumberOfCases),
           value = ILINumberOfCases) %>%
    ungroup() 
  
  df %>% 
    ggplot(aes(x=t_vw, y=value)) + 
    geom_line() + #+ scale_y_log10()
    facet_wrap(~Season, scales ="free_y")
  
  df %>% 
    ggplot(aes(x=t_vw, y=value)) + 
    geom_line() + #+ scale_y_log10()
    facet_wrap(~Season)
  
  df$value[is.na(df$value)] <- 0
  
  inputData <- df %>% 
    filter(!(Season %in% c("2009/2010", "2020/2021", "2021/2022"))) %>%
    group_by(Season) %>%
    mutate(total_inf = sum(value)) %>%
    ungroup() %>% #View()
    filter(total_inf > 0)
  
  inputData %>% 
    ggplot(aes(x=t_vw, y=value)) + 
    geom_line() + #+ scale_y_log10()
    facet_wrap(~Season, scales ="free_y") -> fig0
  print(fig0)
  
  ggsave(paste0("./code/flu/DK_model/Fig_flu_seasons_",country_long,".jpg"), fig0)
  
  # inputData %>%
  #   group_by(Season) %>%
  #   left_join(as_tibble(1:50) %>% rename(t_vw = value), ) %>% 
  #   ungroup %>%
  #   View()
  # 
  # inputData %>% filter(Season == "2012/2013") %>%
  #   right_join(as_tibble(bind_cols(1:50, "2012/2013")) %>% rename(t_vw = ...1, Season = ...2), ) %>% 
  #   View()
  
  
  # inputData %>% 
  #   ggplot(aes(x=t_vw, y=value, color=Season)) + 
  #   geom_line() #+ scale_y_log10()
  
  # inputData$value[is.na(inputData$value)] <- 0
  
  # Load MEM values
  mem_val <- read.csv(paste0("./code/Flu/MEM model/MEMOutput/Baseline_", country_short,".csv"))
  mem_threshold <- mem_val$PreEpidemicThreshold
  if (mem_val$NominatorType != "ILI"){
    stop("MEM threshold is not based on ILI")
  }
  
  # Clean data - use only values after the mem epidemic start
  df_ILI_EU <- inputData %>% 
    group_by(Season) %>%
    mutate(value_infl = value,
           value_cumsum = cumsum(value_infl),
           total_value = max(value_cumsum),
           ind1 = which(value_infl>mem_threshold/2)[1],
           ind2 = which(value_cumsum>0.1*total_value)[1] ) %>%
    filter(t_vw >= ind1) %>%
    mutate(t_vw = 1:length(value_infl)) %>%
    ungroup()
  
  df_ILI_EU %>%
    ggplot(aes(x=t_vw, y=value, color=Season)) + 
    geom_line() 
  
  df_ILI_EU %>%
    ggplot(aes(x=t_vw, y=value)) + 
    geom_line() + 
    facet_wrap(~Season, scales ="free_y")
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Pre-fitting preparation and fit ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  df_param <- NULL
  df_param_FIT <- NULL
  
  for (Y in unique(df_ILI_EU$Season)){
    print(Y)
    if (Y== "2022/2023"){
      next
    }
    
    df_fit <- df_ILI_EU %>%
      filter(Season == Y) 
    
    df_fit %>% ggplot(aes(x=t_vw,y=value_infl)) + geom_point(alpha=0.2) + 
      geom_point(alpha=1)+
      labs(subtitle="Historical wave") 
    
    # #
    # df_new <- df_ILI_EU %>%
    #   filter(Season == "S2015-2016") 
    # #
    
    # ---- |-Setup: cmdstanr ----
    options(mc.cores = detectCores()-1 )
    set_cmdstan_path(path = NULL)
    #mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart.stan') # This compiles the script
    #mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v2.stan') # This compiles the script
    mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v4.stan') # This compiles the script
    
    
    # ---- |-Fit: Historical EU wave  ----
    pop_fit = sum( obtain_demography(country_long) )
    pop_target = pop_fit
    df_fit_0 = df_fit
    wave_start = df_fit
    # Prepare list with parameters for fitting
    stan_list = list(
      n_week_full = nrow(df_fit_0),
      severe_obs_full = (df_fit_0$value_infl) %>% replace_na(replace = 0) %>% round(0),
      n_week_start = nrow(wave_start),
      severe_obs_start = (wave_start$value_infl) %>% replace_na(replace = 0) %>% round(0),
      pop_full=pop_fit,
      pop_start=pop_target
    )
    # Fit
    fit02 <- mod2$sample(
      data = stan_list,
      seed = 12,
      chains = 8,
      parallel_chains = 8,iter_sampling=1500,thin=10,max_treedepth = 15
    )
    
    
    # ---- |-Check the fit results  ----
    
    # create table of parameters
    fit02 %>% gather_draws(SIR_ini[state],
                           I_ini_start[scenario],
                           prop_severe,
                           pop_infect,
                           #prop_severe_start,
                           pop_infect_start[scenario]) %>% 
      mean_qi() -> xp; xp
    
    fit02 %>% gather_draws(SIR_ini[state],
                           prop_severe) -> x
    x$year <- Y
    df_param_FIT <- bind_rows(df_param_FIT, x)
    
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
    round( 1/xp[xp$.variable=="prop_severe",]$.value,1) -> mprob_severe
    (xp[xp$.variable=="SIR_ini"&xp$state==2,".value"]) %>% round(4) -> mI_ini
    (xp[xp$.variable=="SIR_ini"&xp$state==1,".value"]) %>% round(2) -> mS_ini
    (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    fit02 %>% gather_draws(gen_severe_obs_full[t_vw]) %>% 
      mean_qi() %>% left_join(df_fit,by="t_vw") %>% 
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      geom_point(aes(y=value_infl),col="black")+
      labs(subtitle = paste(country_long, Y, ": fit |",
                            "prob_severe:", mprob_severe,"\n",
                            "| S_ini:",mS_ini,
                            "| I_ini:",mI_ini,
                            "| prop inf:",mProp_inf)) -> p_cf0; #p_cf0
    
    #round(xp[xp$.variable=="prop_severe_start"&xp$scenario==1,".value"],3) -> mprob_severe
    (xp[xp$.variable=="I_ini_start"&xp$scenario==1,".value"]) %>% round(4)-> mI_ini
    fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==1) %>% 
      mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
      ggplot(aes(x=t_vw)) + 
      geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
      geom_line(aes(y=.value)) +
      #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
      geom_point(data=wave_start,aes(y=value_infl ),col="darkred") +
      #scale_y_log10() + 
      labs(subtitle = paste("EU scenario: baseline |","\n",
                            #"prob_severe:", mprob_severe,
                            "| I_ini:",mI_ini) ) -> p_cf1; #p_cf1
    mexplain = c("Fitting one EU wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
               prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
               Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
               More transmissible: fit with higher beta
               More vaccination: fit with lower S_ini,
               grey bands showing 80% credible intervals")
    ( ( (p_cf0+p_cf1) ) + plot_annotation(caption = mexplain) )-> ppatch; 
    print(ppatch)
    #ggsave_as(ppatch,"real_aus_eu",30,30)
    
    
    
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ### Save parameters from the historical fit ##########
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    params_hist <- xp
    S_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==1,".value"]
    I_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==2,".value"]
    R_init_hist <- params_hist[params_hist$.variable=="SIR_ini" & params_hist$state==3,".value"]
    prop_severe_hist <- params_hist[params_hist$.variable=="prop_severe",".value"]
    
    # Compare cumulative between data and model
    # Model
    values_model <- fit02 %>% 
      gather_draws(gen_severe_obs[scen,t_vw]) %>% 
      filter(scen==1) %>% 
      mean_qi() %>%
      pull(.value) %>% 
      sum()
    
    # Data
    values_data <- wave_start$value %>% sum()
    
    
    df_tmp <- as_tibble(bind_cols(S_init_hist, I_init_hist, R_init_hist, prop_severe_hist, Y, values_model/values_data, values_model, values_data)) %>% 
      rename(S_init = .value...1, I_init = .value...2, R_init = .value...3, prop_severe = .value...4, season = ...5, prop_total_values = ...6, values_model = ...7, values_data = ...8)
    
    df_param <- bind_rows(df_param, df_tmp)
    
    # print("Going next")
    # next
    # print("------Going next")
    # 
    # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # ### Fit to other seasons ##########
    # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # 
    # Y_new = "S2016-2017"
    # df_fit_new <- df_ILI_EU %>%
    #   filter(Season == Y_new) 
    # 
    # df_fit_new %>% ggplot(aes(x=t_vw,y=value_infl)) + geom_point(alpha=0.2) + 
    #   geom_point(alpha=1)+
    #   labs(subtitle="Historical wave") 
    # 
    # #mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v5.stan') # This compiles the script
    # mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v6.stan') # This compiles the script
    # 
    # 
    # # ---- |-Fit: Historical EU wave  ----
    # pop_fit = sum( obtain_demography(country_long) )
    # pop_target = pop_fit
    # df_fit_0 = df_fit_new
    # wave_start = df_fit_new
    # # Prepare list with parameters for fitting
    # stan_list = list(
    #   n_week_full = nrow(df_fit_0),
    #   severe_obs_full = (df_fit_0$value_infl) %>% replace_na(replace = 0) %>% round(0),
    #   n_week_start = nrow(wave_start),
    #   severe_obs_start = (wave_start$value_infl) %>% replace_na(replace = 0) %>% round(0),
    #   pop_full=pop_fit,
    #   pop_start=pop_target,
    #   S_init_input = df_tmp$S_init,
    #   prop_severe = df_tmp$prop_severe
    # )
    # # Fit
    # fit02 <- mod2$sample(
    #   data = stan_list,
    #   seed = 12,
    #   chains = 8,
    #   parallel_chains = 8,iter_sampling=1500,thin=10,max_treedepth = 15
    # )
    # 
    # 
    # # ---- |-Check the fit results  ----
    # 
    # # create table of parameters
    # fit02 %>% gather_draws(I_init,
    #                        I_ini_start[scenario],
    #                        #prop_severe,
    #                        pop_infect,
    #                        #prop_severe_start,
    #                        pop_infect_start[scenario]) %>% 
    #   mean_qi() -> xp; xp
    # 
    # 
    # fit02 %>% gather_draws(I_ini_logit,
    #                        #prop_severe_logit,
    #                        S_ini_logit,
    #                        I_ini_logit_prior,
    #                        #prop_severe_logit_prior,
    #                        #S_ini_logit_prior
    # ) %>% 
    #   mean_qi() %>% 
    #   ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
    #   geom_pointrange() +
    #   geom_vline(xintercept = c(logit(wave_settings_SIR_Rnull$I_ini),
    #                             logit(wave_settings_SIR_Rnull$prop_severe),
    #                             logit(wave_settings_SIR_Rnull$S_ini)), linetype = 'dotted')
    # 
    # 
    # # ---- |-Plot: Scenario projections  ----
    # #round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
    # round( df_tmp$prop_severe, 3) -> mprob_severe
    # (xp[xp$.variable=="I_init",".value"]) %>% logit() %>% round(1) -> mI_ini
    # round( stan_list$S_init_input, 3) -> mS_ini
    # (xp[xp$.variable=="pop_infect",".value"]) %>% round(2) -> mProp_inf
    # fit02 %>% gather_draws(gen_severe_obs_full[t_vw]) %>% 
    #   mean_qi() %>% left_join(df_fit_0,by="t_vw") %>% 
    #   ggplot(aes(x=t_vw)) + 
    #   geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
    #   geom_line(aes(y=.value)) +
    #   geom_point(aes(y=value_infl),col="black")+
    #   labs(subtitle = paste(country_long, Y_new, ": fit |",
    #                         "prob_severe:", mprob_severe,"\n",
    #                         "| S_ini:",mS_ini,
    #                         "| I_ini:",mI_ini,
    #                         "| prop inf:",mProp_inf)) -> p_cf0; #p_cf0
    # 
    # #round(xp[xp$.variable=="prop_severe_start"&xp$scenario==1,".value"],3) -> mprob_severe
    # (xp[xp$.variable=="I_ini_start"&xp$scenario==1,".value"]) %>% logit() %>% round(1) -> mI_ini
    # fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% filter(scen==1) %>% 
    #   mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
    #   ggplot(aes(x=t_vw)) + 
    #   geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
    #   geom_line(aes(y=.value)) +
    #   #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
    #   geom_point(data=wave_start,aes(y=value_infl ),col="darkred") +
    #   #scale_y_log10() + 
    #   labs(subtitle = paste("EU scenario: baseline |","\n",
    #                         "prob_severe:", mprob_severe,
    #                         "| I_ini:",mI_ini) ) -> p_cf1; #p_cf1
    # mexplain = c("Fitting one EU wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
    #              prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
    #              Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
    #              More transmissible: fit with higher beta
    #              More vaccination: fit with lower S_ini,
    #              grey bands showing 80% credible intervals")
    # ( ( (p_cf0+p_cf1) ) + plot_annotation(caption = mexplain) )-> ppatch; 
    # print(ppatch)
    # #ggsave_as(ppatch,"real_aus_eu",30,30)
    # 
    # 
    # # Compare cumulative between data and model
    # # Model
    # fit02 %>% 
    #   gather_draws(gen_severe_obs[scen,t_vw]) %>% 
    #   filter(scen==1) %>% 
    #   mean_qi() %>%
    #   pull(.value) %>% 
    #   sum()
    # 
    # # Data
    # wave_start$value %>% sum()
    
    
    
    
    
  }
  
  # ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Spectral")
  write_fst(df_param_FIT, paste0("./data/flu_params_",country_long,"_FIT.fst"))
  
  df_param_FIT %>% 
    filter(state==1 | is.na(state)) %>% 
    ungroup() %>%
    select(-state) %>%
    pivot_wider(names_from = .variable, values_from = .value) %>%
    mutate(year = as.factor(year)) %>%
    ggplot(aes(x=SIR_ini, y=prop_severe, color=year)) +
    geom_point() -> fig1
  print(fig1)
  ggsave(paste0("./code/flu/DK_model/Fig_flu_scatter_fit_",country_long,".jpg"), fig1)
  
  
  df_param_FIT %>% 
    filter(state==1 | is.na(state)) %>% 
    ungroup() %>%
    select(-state) %>%
    pivot_wider(names_from = .variable, values_from = .value) %>%
    ggplot(aes(x=SIR_ini, y=prop_severe, group=year, color=year, shape=year)) +
    scale_shape_manual(values=1:nlevels( df_param_FIT$year%>% as.factor() ) ) +
    #labs(title = "Demo more than 6 shapes", x="Theat (deg)", y="Magnitude") +
    #geom_line() + 
    geom_point(size=3)
  
  
  
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ### Fit to the current season and do sensitivity analysis ##########
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ###
  
  ####
  
  # ---- |- Prepare the data  ----
  #Y_new = "S2017-2018"
  df_param_future <- NULL
  #for (Y_new in c("2012/2013", "2015/2016", "2018/2019")){
  for (Y_new in unique(df_ILI_EU$Season)){
    #Y_new = "2015/2016"
    df_fit_new <- df_ILI_EU %>%
      filter(Season == Y_new) #%>% 
      # filter(t_vw>5) %>%
      # group_by(Season) %>%
      # mutate(t_vw = 1:length(ILINumberOfCases)) %>%
      # ungroup() 
    
    df_fit_new %>% ggplot(aes(x=t_vw,y=value_infl)) + geom_point(alpha=0.2) + 
      geom_point(alpha=1)+
      labs(subtitle="Historical wave") 
    
    
    if (Y_new == "2022/2023"){
      df_fit_new_short <- df_fit_new 
    } else {
      df_fit_new_short <- df_fit_new %>% 
        mutate(value_cumsum = cumsum(value),
               total_value = max(value_cumsum)) %>% 
        filter(value_cumsum<0.2*total_value) %>% 
        select(-total_value, -value_cumsum, -ind1, -ind2)
    }
    
    
    df_fit_new_short %>% ggplot(aes(x=t_vw,y=value_infl)) + geom_point(alpha=0.2) + 
      geom_point(alpha=1)+
      labs(subtitle="Historical wave") 
    
    
    
    
    #for (j in 1:nrow(df_param)){
    for (j in 1:10){
      
      # ---- |- Set the parameters  ----
      S_mean = df_param$S_init %>% mean()
      S_sd = df_param$S_init %>% sd()
      S_start0 = rnorm(1, S_mean, S_sd)
      
      prop_severe_mean = df_param$prop_severe %>% mean()
      prop_severe_sd = df_param$prop_severe %>% sd()
      prop_severe_start0 = rnorm(1, prop_severe_mean, prop_severe_sd)
      
      df_param_FIT %>% 
        filter(state==1 | is.na(state)) %>% 
        ungroup() %>%
        select(-state) %>%
        pivot_wider(names_from = .variable, values_from = .value) -> X
      
      k = sample(1:nrow(X),1)
      S_start0 = X$SIR_ini[k]
      prop_severe_start0 = X$prop_severe[k]
      year_sample = X$year[k]
      
      # S_start0 = df_param$S_init[j]
      # prop_severe_start0 = df_param$prop_severe[j]
      
      # Define the vaccination interventions
      n_interventions = 3 # Including baseline!
      interventions_impact = c(0, 0.05, 0.10) # Decrease in S due to vaccination
      vacc_timing = nrow(df_fit_new_short)*7 + 1 # Time of vaccination in days (all vaccine doses given on this day)
      
      # ---- |- Load the model  ----
      mod2 <- cmdstan_model(stan_file='./stan/immu_dyn_02_eustart_v7.stan') # This compiles the script
      
      
      
      # ---- |-Fit  ----
      pop_fit = sum( obtain_demography(country_long) )
      pop_target = pop_fit
      df_fit_0 = df_fit_new
      wave_start = df_fit_new_short
      n_weeks0 = ifelse(Y_new == "2022/2023", 40, nrow(df_fit_0))
      severe_vec = rep(0, n_weeks0)
      # Prepare list with parameters for fitting
      stan_list = list(
        #n_week_full = nrow(df_fit_0),
        n_week_full = n_weeks0,
        #severe_obs_full = 0*(df_fit_0$value_infl) %>% replace_na(replace = 0) %>% round(0),
        severe_obs_full = severe_vec,
        n_week_start = nrow(wave_start),
        severe_obs_start = (wave_start$value_infl) %>% replace_na(replace = 0) %>% round(0),
        pop_full=pop_fit,
        pop_start=pop_target,
        S_init_input = S_start0,
        prop_severe = prop_severe_start0,
        n_interventions = n_interventions,
        interventions_impact = interventions_impact,
        vacc_timing = vacc_timing
      )
      # Fit
      fit02 <- mod2$sample(
        data = stan_list,
        seed = 12,
        chains = 8,
        parallel_chains = 8,iter_sampling=1500,thin=10,max_treedepth = 15
      )
      
      
      # ---- |-Check the fit results  ----
      
      # create table of parameters
      fit02 %>% gather_draws(#I_init,
        # I_ini_start[scenario],
        I_ini_start,
        prop_severe_start,
        #prop_severe,
        pop_infect,
        #prop_severe_start,
        pop_infect_start[scenario]) %>% 
        mean_qi() -> xp; xp
      
      
      fit02 %>% gather_draws(I_ini_logit,
                             #prop_severe_logit,
                             S_ini_logit,
                             I_ini_logit_prior,
                             #prop_severe_logit_prior,
                             #S_ini_logit_prior
      ) %>% 
        mean_qi() %>% 
        ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper)) +
        geom_pointrange() +
        geom_vline(xintercept = c(logit(wave_settings_SIR_Rnull$I_ini),
                                  logit(wave_settings_SIR_Rnull$prop_severe),
                                  logit(wave_settings_SIR_Rnull$S_ini)), linetype = 'dotted')
      
      
      # Obtain cumulative burden for model (and compare to data)
      # Model
      L <- NULL
      for (jj in 1:n_interventions){
        x <- fit02 %>% 
          gather_draws(gen_severe_obs[scen,t_vw]) %>% 
          filter(scen==jj) %>% 
          mean_qi() %>%
          pull(.value) %>% 
          sum()
        L <- append(L, x)
        
      }
      values_total_model <- L
      
      # Data
      values_total_data <- df_fit_0$value %>% sum()
      
      
      # ---- |-Plot: Scenario projections  ----
      if (T){
        #round(xp[xp$.variable=="prop_severe",".value"],3) -> mprob_severe
        round( prop_severe_start0, 5) -> mprob_severe
        (xp[xp$.variable=="I_init",".value"]) %>% logit() %>% round(1) -> mI_ini
        round( stan_list$S_init_input, 3) -> mS_ini
        #(xp[xp$.variable=="prop_severe",".value"]) %>% round(2) -> mProp_inf
        fit02 %>% gather_draws(gen_severe_obs_full[t_vw]) %>% 
          mean_qi() %>% left_join(df_fit_0,by="t_vw") %>% 
          ggplot(aes(x=t_vw)) + 
          geom_ribbon(aes(ymin=.lower,ymax=.upper)) + 
          geom_line(aes(y=.value)) +
          geom_point(aes(y=value_infl),col="black")+
          labs(subtitle = paste(country_long, Y_new, ": fit |",
                                "prob_severe:", mprob_severe,"\n",
                                "| S_ini:",mS_ini,
                                "| prop inf:",mProp_inf,
                                "| I_ini:",mI_ini)) -> p_cf0; #p_cf0
        
        #round(xp[xp$.variable=="prop_severe_start"&xp$scenario==1,".value"],3) -> mprob_severe
        (xp[xp$.variable=="I_ini_start",".value"]) %>% logit() %>% round(1) -> mI_ini
        #(xp[xp$.variable=="pop_infect_start",".value"]) %>% round(2) -> mProp_inf
        ( values_total_model/(xp[xp$.variable=="prop_severe_start",]$.value) / pop_fit) %>% round(2) -> mProp_inf
        fit02 %>% gather_draws(gen_severe_obs[scen,t_vw]) %>% #filter(scen==1) %>% 
          mean_qi() %>% left_join(wave_start,by="t_vw") %>% 
          ggplot(aes(x=t_vw)) + 
          geom_ribbon(aes(ymin=.lower,ymax=.upper, color=as.factor(scen))) + 
          geom_line(aes(y=.value, color=as.factor(scen))) +
          #geom_point(aes(y=severe_obs_weekly),col="lightgrey")+
          geom_point(data=df_fit_0,aes(y=value_infl )) +
          geom_point(data=wave_start,aes(y=value_infl ),col="darkred") +
          #scale_y_log10() + 
          labs(subtitle = paste("EU scenario: baseline |","\n",
                                "prob_severe:", mprob_severe,
                                "| I_ini:",mI_ini,
                                "| prop inf:",mProp_inf) ) -> p_cf1; #p_cf1
        
        
        
        mexplain = c("Fitting one EU wave to SIR model with parameters I_ini (initially infected), S_ini (initially susceptible), 
               prob_severe (proportion of infected observed, dark factor relative to fully-immunising infections), 
               Baseline: Fitting EU wave with own I_ini (different wave timing) and prob_severe (different surveillance), but the same S_ini. 
               More transmissible: fit with higher beta
               More vaccination: fit with lower S_ini,
               grey bands showing 80% credible intervals")
        ( ( (p_cf0+p_cf1) ) + plot_annotation(caption = mexplain) )-> ppatch; 
        print(ppatch)
        #ggsave_as(ppatch,"real_aus_eu",30,30)
      }
      
      # ---- |- Save the relevant parameters  ----
      
      # col = bind_cols(j, S_start0, prop_severe_start0, mProp_inf, Y_new, year_sample)
      # for (jj in 1:n_interventions){
      #   x1 <- values_total_model[jj] / values_total_data 
      #   col <- bind_cols(col, x1)
      # }
      # for (jj in 1:n_interventions){
      #   x2 <- values_total_model[jj] / values_total_model[1]
      #   col <- bind_cols(col, x2)
      # }
      
      ( values_total_model/(xp[xp$.variable=="prop_severe_start",]$.value) / pop_fit) %>% round(2) -> mProp_inf
      
      df_tmp <- as_tibble(bind_cols(j, S_start0, prop_severe_start0, mProp_inf, Y_new, values_total_model/values_total_model[1], year_sample, values_total_model/values_total_data, (1:n_interventions), country_long, values_total_data )) %>% 
        rename( iteration = ...1, S_init = ...2, prop_severe = ...3, prop_inf = ...4, Season = ...5, intervention_impact = ...6, year_sample = ...7, prop_total_values = ...8, scen = ...9, Country = ...10, values_total_data = ...11  )
      
      df_param_future <- bind_rows(df_param_future, df_tmp)
      
      
    }
  }
  
  write_fst(df_param_future, paste0("./data/flu_params_",country_long,"_future.fst"))
  
  
}





# Load the data
X1 <- read_fst(paste0("./data/flu_params_","Latvia","_future.fst"))
X2 <- read_fst(paste0("./data/flu_params_","Slovenia","_future.fst"))
X3 <- read_fst(paste0("./data/flu_params_","Austria","_future.fst"))

df_param_future <- bind_rows(X1,X2,X3)



# Plot S_init and prop_severe vs true burden (from data)

df_param_future %>% 
  group_by(Country,Season,scen, values_total_data) %>%
  summarise(m = mean(S_init ),
            sd = sd(S_init)) %>%
  ungroup() %>%
  ggplot(aes(x=m, y=values_total_data)) + 
  geom_point() + 
  scale_colour_gradient(low="white", high = "black") +
  facet_wrap(~Country, scales = "free")


df_param_future %>% 
  group_by(Country,Season,scen, values_total_data) %>%
  summarise(m = mean(prop_severe),
            sd = sd(prop_severe)) %>%
  ungroup() %>%
  ggplot(aes(x=m, y=values_total_data)) + 
  geom_point() + 
  scale_colour_gradient(low="white", high = "black") +
  facet_wrap(~Country, scales = "free")



# Intervention impact
df_param_future %>% filter(Country=="Slovenia") %>%
  group_by(Country, scen, Season) %>%
  mutate(burd = ifelse(values_total_data<400, 1, 2)) %>%
  ungroup() %>%
  mutate(Scenario = as.factor(scen)) %>%
  mutate(Scenario = ifelse(Scenario==1, "Scenario 0 - Baseline",
                           ifelse(Scenario==2, "Scenario 1 - 5% decrease in S", 
                                  "Scenario 2 - 10% decrease in S"))) %>%
  ggplot(aes(x=Season, y=1-intervention_impact, fill=as.factor(Scenario))) + 
  labs(y="Relative decrease in seasonal burden") +
  geom_boxplot() #+ 
  #facet_wrap(~Scenario)

df_param_future %>% #filter(Country=="Slovenia") %>%
  mutate(Scenario = as.factor(scen)) %>%
  mutate(Scenario = ifelse(Scenario==1, "Scenario 0 - Baseline",
                           ifelse(Scenario==2, "Scenario 1 - 5% decrease in S", 
                                  "Scenario 2 - 10% decrease in S"))) %>%
  ggplot(aes(x=Country, y=1-intervention_impact, fill=as.factor(Scenario))) + 
  labs(y="Relative decrease in seasonal burden") +
  geom_boxplot() 




#
df_param_future %>% filter(Season == "2022/2023") %>% head()



#
df_param_future %>% 
  ggplot(aes(x=prop_severe, y=intervention_impact)) + 
  geom_point(aes(fill=S_init)) + 
  scale_colour_gradient(low="white", high = "black")



df_param_future %>% 
  ggplot(aes(x=S_init, y=intervention_impact, color=year_sample)) + 
  geom_point(aes(size=prop_severe)) +
  facet_wrap(~scen, scales ="free_y")

df_param_future %>% 
  ggplot(aes(x=prop_severe, y=intervention_impact, color=year_sample)) + 
  geom_point(aes(size=S_init)) +
  facet_wrap(~scen, scales ="free_y")

df_param_future %>% 
  ggplot(aes(x=S_init, y=prop_total_values, color=year_sample)) + 
  geom_point(aes(size=prop_severe)) +
  facet_wrap(~scen, scales ="free_y")

df_param_future %>% 
  ggplot(aes(x=prop_severe, y=prop_total_values, color=year_sample)) + 
  geom_point(aes(size=S_init)) +
  facet_wrap(~scen, scales ="free_y")

df_param_future %>% 
  ggplot(aes(x=1-intervention_impact, color=as.factor(scen))) + 
  geom_histogram() 

df_param_future %>% 
  ggplot(aes(x=1-intervention_impact, color=as.factor(year_sample))) + 
  geom_histogram() +
  facet_wrap(~scen, scales ="free_y")



df_param_future %>% 
  ggplot(aes(x=prop_inf, color=year_sample)) + 
  geom_histogram() +
  facet_wrap(~scen, scales ="free_y")

df_param_future %>% 
  ggplot(aes(x=prop_total_values)) + 
  geom_histogram() +
  facet_wrap(~scen, scales ="free_y")


#
df_param_future %>% 
  ggplot(aes(y = values_total_data, x=1-intervention_impact)) +
  geom_point()

df_param_future %>%
  filter(scen != 1) %>% 
  group_by(Season, scen, values_total_data, Country) %>%
  summarise(m = mean(1-intervention_impact),
            s = sd(1-intervention_impact)) %>%
  ungroup() %>%
  ggplot(aes(x=m, y=values_total_data, color=as.factor(scen))) +
  geom_point() +
  geom_errorbar(aes(xmin=m-s, xmax=m+s), width=.2,
                position=position_dodge(.9)) + 
  facet_wrap(Country~scen, scales ="free")


df_param_future %>%
  filter(scen != 1) %>% 
  group_by(Season, scen, values_total_data, Country) %>%
  summarise(m = mean(1-intervention_impact),
            s = sd(1-intervention_impact)) %>%
  ungroup() %>%
  ggplot(aes(x=m, y=values_total_data, color=as.factor(scen))) +
  geom_point() + 
  facet_wrap(Country ~ scen, scales ="free")


# BOXPLOT
df_param_future %>% 
  ggplot(aes(x=Season, y=1-intervention_impact, fill=as.factor(scen))) + 
  geom_boxplot() #+
# facet_wrap(~scen)


df_param_future %>% 
  ggplot(aes(x=scen, y=1-intervention_impact)) + 
  geom_violin() +
  facet_wrap(~scen)

