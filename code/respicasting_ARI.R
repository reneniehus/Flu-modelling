# ---- |-Set up ----
source("code/setup.R")
# ---- |-Settings and parameters ----
clc()
indicator = "ARI"
myorigin = c("2023-12-27") # last: c("2023-12-20")
truth_date_latest = ymd("2023-12-17") # last: ymd("2023-12-10")
mytarget = paste0(indicator," incidence")
weeks_forecast = 6
data_points_fitted = 8
size_group = 2
myeps = 1/100000
myquantiles = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
fit_rerun = T

tranv_v = c("log","sqrt")
for (i_trans in tranv_v) { # run same model for diff scaling, i_trans=tranv_v[1]
  mytransformation = i_trans
  
  # Transformed settings and parameters
  myorigin_date = ymd(myorigin)
  x_linear_adjust = sum( seq(1:size_group)-1 ) / size_group
  if (mytransformation == "log") {
    data_transform = base::log
    data_backtransform = base::exp
    cat("Transformation: log")
  }
  if (mytransformation == "sqrt") {
    data_transform = base::sqrt
    data_backtransform = function(x){x^2}
    cat("Transformation: sqrt")
  }
  iter_stan = 2000
  
  # ---- |-Load data for all countries ----
  # check data
  # https://github.com/european-modelling-hubs/ari-forecast-hub/tree/main/target-data/ERVISS
  data2 = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/ari-forecast-hub/main/target-data/ERVISS/latest-ARI_incidence.csv",show_col_types=F)
  data = data2
  #browser()
  # ---- |-Potentially remove countries ----
  countries_select = data %>% group_by(location) %>% 
    summarise(max_truth_date = max(truth_date)) %>% 
    filter(max_truth_date %in% c(truth_date_latest,truth_date_latest-7) ) %>% pull(location)
  countries_removed = data %>% select(location) %>% distinct() %>% filter(!location %in% countries_select) %>% nrow()
  data = data %>%  filter(location %in% countries_select)
  country_v = unique(data$location) ;
  pr=paste("Data for",countries_removed,"locations not up to date\n"); cat(red(pr))
  pr=paste("Data loaded for",length(country_v),"locations\n"); cat(blue(pr))
  
  # ---- |-Loop: for each country ----
  respicast_df_list = list()
  for (country_i in country_v) {
    # loop helpers
    population_i = aux$pops %>% filter(location==country_i) %>% pull(population)
    
    # ---- |-Filtering  ----
    sentence=paste0(">Running: ",EU_long(country_i)," (pop:",population_i,")\n");cat(sentence)
    data_loc = data %>% 
      filter(location==country_i) %>% 
      arrange(truth_date) %>% 
      mutate(day_diff=as.numeric(truth_date-lag(truth_date,1)),.before = location) %>% 
      tail(n=data_points_fitted) # filter only the last data_points_fitted weekly data point
    if( any(data_loc$day_diff[-1]!=7) ) { pr=paste("Warning: Check day diffs for:",country_i,"\n"); cat(red(pr))}
    if (F) data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point() 
    
    # ---- |-Data transform  ----
    data_loc = data_loc %>% 
      # 1: first transform to ensure not exceeding pop size -> to odds
      mutate(value_odds = odds(value/population_i) ) %>% 
      # 2: model-specific tranform
      mutate(value_transformed = data_transform( value_odds %>% zero_plus_eps( eps=myeps ) )) 
    
    # ---- |-Data&other warnings  ----
    if ( any(is.na(data_loc$value))  ) {
      pr=paste("Warning: NA in values:",country_i,"\n"); cat(red(pr))
      next; # don't run rest of for-loop, go to next country
    }
    
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
    ) %>% 
      mutate(x_linear=x_linear-x_linear_adjust )
    
    if (F) {
      data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
      data_loc %>% ggplot(aes(x=truth_date,y=value_transformed)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
    }
    # filter out points without full group
    df_stan = data_loc %>% filter(belongs_complete_group==1)
    
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
      df_stan_predict_extra[[wahead]]$value <- df_stan_predict_extra[[wahead]]$value_transformed <- NA
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
      
      group_intercept= df_stan %>% group_by(point_group) %>% summarise( group_intercept=mean(value_transformed) ) %>% pull(group_intercept),
      x_linear = df_stan_predict_all$x_linear,
      y = df_stan$value_transformed,
      #
      prior_intercept_sd=0.1,
      #prior_slope_sd=1,
      prior_slope_diff=0.05
    )
    # ---- |-Fit stan model ----
    if (mytransformation=="sqrt") myfile=paste0("../Big data/respicasting_",country_i,"sqrt_fit01",indicator,".Rdata" )
    if (mytransformation=="log") myfile=paste0("../Big data/respicasting_",country_i,"log_fit01",indicator,".Rdata" )
    
    if (fit_rerun==T) {
      mod1_path = c("./stan/piecewise_01_starting.stan")
      options(mc.cores = 8 )
      fit01=rstan::stan(
        file=mod1_path,
        chains=8 ,thin=8,iter=iter_stan,
        seed=12, cores = getOption("mc.cores", 1L),
        open_progress = F,
        control=list(
          #adapt_delta=0.9,
          max_treedepth=14
        ),
        data=stan_list
      ) # 
      save(fit01,stan_list,file=myfile )
    }
    load(file=myfile )
    fit = fit01
    # ---- |-Checks by hand ----
    if (F) { # checks by hand
      
      fit@date # Dec  4 21:26:58 2023
      fit@model_pars
      
      mypars = c("slope","slope_diff")
      precis(fit,mypars,depth=2)
    }
    
    if (F) { # add predictions and look
      data_mod = fit %>% gather_draws( gen_y[df_i] ) %>% 
        mutate(.value=data_backtransform(.value)) %>% 
        mutate(.value=inv_odds(.value)*population_i) %>%
        mode_qi(.width=0.5) %>% select(df_i ,.value,.lower,.upper)
      
      data_res = df_stan_predict_all %>% left_join(data_mod , by="df_i")
      
      data_res %>% ggplot(aes(x=truth_date)) + 
        geom_point(aes(y=value)) +
        geom_line(aes(y=.value)) + 
        geom_ribbon(aes(ymin=.lower, ymax=.upper ))
    }
    
    # ---- |-format data for submission ----
    # https://github.com/european-modelling-hubs/flu-forecast-hub/wiki/Submission-format
    df_n_data = nrow(df_stan)
    
    # only point estimate
    model_gen01 = fit %>% gather_draws( gen_y[df_i] ) %>% 
      filter(df_i > df_n_data ) %>% 
      mutate(.value=data_backtransform(.value)) %>% 
      mutate(.value=inv_odds(.value)*population_i) %>%
      mode_qi(.width=c(0.5) ) %>% 
      select(df_i ,.value) %>% 
      mutate(output_type="median",output_type_id="") %>% 
      select(df_i,value=.value,output_type,output_type_id); # g(model_gen01)
    # quanitles
    model_gen02 = fit %>% gather_draws( gen_y[df_i] ) %>% 
      filter(df_i > df_n_data ) %>% 
      mutate(.value=data_backtransform(.value)) %>% 
      mutate(.value=inv_odds(.value)*population_i) %>%
      column_stats_ingroups(mycolumn=.value,mygroup = "df_i",probs=myquantiles ) %>% 
      mutate(output_type="quantile",output_type_id=as.character(round(quant,3))) %>% 
      select(df_i,value=val,output_type,output_type_id); # g(model_gen02)
    # put it all together
    respicast_df_list[[country_i]] = bind_rows(model_gen01 , model_gen02) %>% left_join( df_stan_predict_all %>% rename(value_truth=value) , by="df_i" ) %>% 
      mutate(
        origin_date=as.character(myorigin_date),
        target=mytarget,
        horizon= as.numeric((truth_date-(myorigin_date-3))/7 + 1) ,
        truth_date = as.character(truth_date)
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
  
  # filter correct forecasting dates
  respicast_df_unfiltered = respicast_df
  respicast_df = respicast_df %>% filter_log(horizon%in%c(1:4)) 
  
  
  if (mytransformation=="log") write_csv( respicast_df , file=paste0(
    "./output/ari-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_green",indicator,".csv") )
  if (mytransformation=="sqrt") write_csv( respicast_df , file=paste0(
    "./output/ari-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_blue",indicator,".csv") )
  
} # through different scaling

respicast_df %>% select(horizon,target_end_date) %>% distinct()

# https://github.com/Infectious-Disease-Modeling-Hubs/hubVis
# remotes::install_github("Infectious-Disease-Modeling-Hubs/hubVis")
if (F){
  
  library(hubVis)
  mod_log = read_csv(
    file=paste0(
      "./output/ari-forecast-hub/ECDC-norrsken_green/",myorigin,"-ECDC-norrsken_green",indicator,".csv"),show_col_types=F )
  mod_sqrt = read_csv(
    file=paste0(
      "./output/ari-forecast-hub/ECDC-norrsken_blue/",myorigin,"-ECDC-norrsken_blue",indicator,".csv"),show_col_types=F )
  plot_mod_log = mod_log %>% mutate(model_id="log",
                                    target_date=target_end_date,
                                    output_type_id = as.numeric(output_type_id)) %>% 
    filter(output_type != "median")
  plot_mod_sqrt = mod_sqrt %>% mutate(model_id="sqrt",
                                      target_date=target_end_date,
                                      output_type_id = as.numeric(output_type_id)) %>% 
    filter(output_type != "median")
  
  
  plot_step_ahead_model_output(bind_rows(plot_mod_log,plot_mod_sqrt),
                               data %>% mutate(time_idx=truth_date) %>% filter(truth_date>ymd("2023-11-15")),
                               facet=c("location"), facet_scales = "free",
                               intervals = c(0.5),interactive=F) 
  
}



