# load the fake date
data = read_csv(file="./data/2020-03-08-ILI_incidence.csv",show_col_types=F)
data = data %>% 
  mutate(value_log = log(value)) %>% 
  mutate(value_sqrt = sqrt(value))

# filter country
data_loc = data %>% filter_log(location=="IT") 
data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point() 

# make it more tricky
data_loc = data_loc %>% slice(2:10)

# add grouping
point_n = nrow(data_loc)
group_n = round(point_n/2)
group_left = point_n%%2

point_identities_complete_group = rep(1:group_n,each=2)
point_identities_incomplete_group = rep(0,each=group_left)
point_identities = c( point_identities_incomplete_group , point_identities_complete_group  )
point_identities = point_identities + 1

x_linear = c( 
  seq_along(point_identities_incomplete_group) , 
  rep(1:2,group_n)  )
x_linear = x_linear-1

belongs_complete_group = c( rep(0,each=group_left) , rep(1,each = length(point_identities_complete_group)  )  )

data_loc = data_loc %>% mutate( 
  point_group = point_identities, 
  belongs_complete_group=belongs_complete_group,  
  x_linear = x_linear
  )

data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
data_loc %>% ggplot(aes(x=truth_date,y=value_log)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 
data_loc %>% ggplot(aes(x=truth_date,y=value_sqrt)) + geom_line() + geom_point(aes(col=as.factor(point_group),pch=as.factor(belongs_complete_group)),size=7) 

# filter
data_loc_for_stan = data_loc %>% filter_log(belongs_complete_group==1)

# extend stan_df for predictions
data_loc_pred = data_loc_for_stan
data_loc_pred_last = data_loc_pred[nrow(data_loc_pred),]
# 1 week pred
wahead = 1
data_loc_pred_1w = data_loc_pred_last
# same: location, issue_date
data_loc_pred_1w$truth_date = data_loc_pred_last$truth_date + 7
data_loc_pred_1w$value <- data_loc_pred_1w$value_log <- data_loc_pred_1w$value_sqrt <- NA
data_loc_pred_1w$point_group = data_loc_pred_last$point_group + 1
data_loc_pred_1w$x_linear = wahead

data_loc_pred_last = data_loc_pred_1w
# 2 week pred
wahead = 2
data_loc_pred_2w = data_loc_pred_last
data_loc_pred_2w$truth_date = data_loc_pred_last$truth_date + 7
data_loc_pred_2w$value <- data_loc_pred_2w$value_log <- data_loc_pred_2w$value_sqrt <- NA
data_loc_pred_2w$point_group = data_loc_pred_last$point_group + 0
data_loc_pred_2w$x_linear = wahead

data_loc_pred_last = data_loc_pred_2w

# 3 week pred
wahead = 3
data_loc_pred_3w = data_loc_pred_last
data_loc_pred_3w$truth_date = data_loc_pred_last$truth_date + 7
data_loc_pred_3w$value <- data_loc_pred_3w$value_log <- data_loc_pred_3w$value_sqrt <- NA
data_loc_pred_3w$point_group = data_loc_pred_last$point_group + 0
data_loc_pred_3w$x_linear = wahead

data_loc_pred_last = data_loc_pred_3w

# 4 week pred
wahead = 4
data_loc_pred_4w = data_loc_pred_last
data_loc_pred_4w$truth_date = data_loc_pred_last$truth_date + 7
data_loc_pred_4w$value <- data_loc_pred_4w$value_log <- data_loc_pred_4w$value_sqrt <- NA
data_loc_pred_4w$point_group = data_loc_pred_last$point_group + 0
data_loc_pred_4w$x_linear = wahead

data_loc_pred_last = data_loc_pred_4w

data_loc_pred_all = bind_rows(
  data_loc_pred,
  data_loc_pred_1w,data_loc_pred_2w,data_loc_pred_3w,data_loc_pred_4w
) 


# make stan list
stan_list = list(
  n = nrow(data_loc_for_stan),
  n_predict = 4,
  n_with_predict = nrow(data_loc_pred_all),
  n_group = n_distinct(data_loc_for_stan$point_group),
  group = fct_inorder(as.character(data_loc_for_stan$point_group)) %>% as.numeric(),
  
  group_intercept=data_loc_for_stan %>% filter(x_linear==0) %>% pull(value_log),
  x_linear = data_loc_pred_all$x_linear,
  y = data_loc_for_stan$value_log,
  #
  prior_intercept_sd=0.1,
  prior_slope_sd=0.1,
  prior_slope_diff_sd=0.1
)

# fit stan model
mod1_path = c("./stan/piecewise_01_starting.stan")
options(mc.cores = 8 )
fit01=rstan::stan(
  file=mod1_path,
  chains=8 ,thin=8,iter=150,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    max_treedepth=14
  ),
  data=stan_list
) # 
save(fit01,stan_list,file=paste0("../Big data/respicasting_fit01.Rdata") )
load(file=paste0("../Big data/respicasting_fit01.Rdata") )

# checks by hand
fit = fit01
fit@date
fit@model_pars

mypars = c("mu")
precis(fit,mypars,depth=2)

# adding predictions
data_loc_pred_all

