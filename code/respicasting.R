# ---- |-Set up ----
source("code/setup.R")
# ---- |-Stuff for settings ----
data_points_fitted = 10
size_group = 2
weeks_forecast = 4
myorigin_date = ymd("2020-03-08")
mytarget = "ILI incidence"
myquantiles = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
fit_rerun = T

# ---- |-Load data for all countries ----
data = read_csv(file="./data/2020-03-08-ILI_incidence.csv",show_col_types=F)
data = data %>% 
  mutate(value_log = log( value %>% zero_plus_eps( eps=1/100000 ) )) %>% 
  mutate(value_sqrt = sqrt( value %>% zero_plus_eps( eps=1/100000 )  ))

# ---- |-Loop: for each country ----
country_v = unique(data$location) ; 
respicast_df_list = list()
for (country_i in country_v) {
  # ---- |-Filtering  ----
  sentence=paste0("Running: ",country_i);print(sentence)
  data_loc = data %>% 
    filter_log(location==country_i) %>% 
    tail(n=data_points_fitted) # filter only the last data_points_fitted weekly data point
  if (F) data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point() 
  
  # ---- |-Adding grouping to data points  ----
  point_n = nrow(data_loc)
  group_n = round(point_n/size_group)
  group_left = point_n%%size_group
  
  point_identities_complete_group = rep(1:group_n,each=size_group)
  point_identities_incomplete_group = rep(0,each=group_left)
  point_identities = c( point_identities_incomplete_group , point_identities_complete_group  )
  point_identities = point_identities + 1
  
  x_linear = c( 
    seq_along(point_identities_incomplete_group) , 
    rep(1:size_group,group_n)  )
  x_linear = x_linear-1
  
  belongs_complete_group = c( rep(0,each=group_left) , rep(1,each = length(point_identities_complete_group)  )  )
  
  data_loc = data_loc %>% mutate( 
    point_group = point_identities, 
    belongs_complete_group=belongs_complete_group,  
    x_linear = x_linear
  )
  if (F) {
  data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
  data_loc %>% ggplot(aes(x=truth_date,y=value_log)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
  data_loc %>% ggplot(aes(x=truth_date,y=value_sqrt)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
  }
  # filter out points without full group
  df_stan = data_loc %>% filter_log(belongs_complete_group==1)
  
  # ---- |-Extend data frame with "future rows" ----
  #
  df_stan_predict = df_stan
  df_stan_predict_last = df_stan_predict[nrow(df_stan_predict),]
  point_group_predict=df_stan_predict_last$point_group + 1
  df_stan_predict_extra = list()
  for (week_ahead_i in 1:weeks_forecast) {
    wahead = week_ahead_i
    df_stan_predict_extra[[wahead]] = df_stan_predict_last
    # stays same: location, issue_date
    df_stan_predict_extra[[wahead]]$truth_date = df_stan_predict_last$truth_date + 7
    df_stan_predict_extra[[wahead]]$value <- df_stan_predict_extra[[wahead]]$value_log <- df_stan_predict_extra[[wahead]]$value_sqrt <- NA
    df_stan_predict_extra[[wahead]]$point_group = point_group_predict
    df_stan_predict_extra[[wahead]]$x_linear = wahead
    
    df_stan_predict_last = df_stan_predict_extra[[wahead]]
  }
  df_stan_predict_extra = bind_rows(df_stan_predict_extra)
  df_stan_predict_all = bind_rows(df_stan_predict,df_stan_predict_extra) %>% 
    mutate(df_i=1:nrow(.) ) # add row counter to merge predictions
  
  # ---- |-Prepare list for stan ----
  # make stan list
  stan_list = list(
    n = nrow(df_stan),
    n_predict = weeks_forecast,
    n_with_predict = nrow(df_stan_predict_all),
    n_group = n_distinct(df_stan$point_group),
    group = fct_inorder(as.character(df_stan$point_group)) %>% as.numeric(),
    
    group_intercept=df_stan %>% filter(x_linear==0) %>% pull(value_log),
    x_linear = df_stan_predict_all$x_linear,
    y = df_stan$value_log,
    #
    prior_intercept_sd=0.1,
    prior_slope_sd=1,
    prior_slope_diff_sd=0.1
  )
  # ---- |-Fit stan model ----
  myfile=paste0("../Big data/respicasting_",country_i,"_fit01.Rdata" )
  if (fit_rerun==T) {
    mod1_path = c("./stan/piecewise_01_starting.stan")
    options(mc.cores = 8 )
    fit01=rstan::stan(
      file=mod1_path,
      chains=8 ,thin=8,iter=1000,
      seed=12, cores = getOption("mc.cores", 1L),
      control=list(
        #adapt_delta=0.9,
        max_treedepth=14
      ),
      data=stan_list
    ) # 
    save(fit01,stan_list,file=myfile )
  }
  load(file=myfile )
  
  # ---- |-Checks by hand ----
  if (F) { # checks by hand
    
    fit = fit01
    fit@date # Nov 12 17:26:41
    fit@model_pars
    
    mypars = c("slope")
    precis(fit,mypars,depth=2)
  }
  
  if (F) { # add predictions and look
    data_mod = fit %>% gather_draws( gen_y[df_i] ) %>% 
      mode_qi(.width=0.5) %>% select(df_i ,.value,.lower,.upper)
    
    data_res = df_stan_predict_all %>% left_join(data_mod , by="df_i")
    
    data_res %>% ggplot(aes(x=truth_date)) + 
      geom_point(aes(y=value_log)) +
      geom_line(aes(y=.value)) + 
      geom_ribbon(aes(ymin=.lower, ymax=.upper ))
  }
  
  # ---- |-format data for submission ----
  # https://github.com/european-modelling-hubs/flu-forecast-hub/wiki/Submission-format
  # Please use 2020-03-08 as origin_date for this test
  df_n_data = nrow(df_stan)
  
  # only point estimate
  model_gen01 = fit %>% gather_draws( gen_y[df_i] ) %>% 
    filter(df_i > df_n_data ) %>% 
    mode_qi(.width=c(0.5) ) %>% 
    select(df_i ,.value) %>% 
    mutate(output_type="median",output_type_id="") %>% 
    select(df_i,value=.value,output_type,output_type_id); g(model_gen01)
  # quanitles
  model_gen02 = fit %>% gather_draws( gen_y[df_i] ) %>% 
    filter(df_i > df_n_data ) %>% 
    column_stats_ingroups(mycolumn=.value,mygroup = "df_i",probs=myquantiles ) %>% 
    mutate(output_type="quantile",output_type_id=format(quant,nsmall=2)) %>% 
    select(df_i,value=val,output_type,output_type_id); g(model_gen02)
  
  respicast_df_list[[country_i]] = bind_rows(model_gen01 , model_gen02) %>% left_join( df_stan_predict_all %>% rename(value_truth=value) , by="df_i" ) %>% 
    mutate(
      origin_date=myorigin_date,
      target=mytarget,
      horizon=df_i-df_n_data,
    ) %>% 
    select( 
      origin_date,
      target,
      target_end_date=truth_date,
      horizon,
      location,
      output_type,
      output_type_id,
      value
    ) 
} # end of country loop
respicast_df = bind_rows(respicast_df_list) 

length(country_v) * (length(myquantiles) + 1 ) * weeks_forecast

write_csv( respicast_df ,file="./output/2020-03-08-ECDC-norrsken.csv" )
