# load the fake date
data = read_csv(file="./data/2020-03-08-ILI_incidence.csv",show_col_types=F)
data = data %>% mutate(value_log = log(value))

# filter country
data_loc = data %>% filter_log(location=="IT") 
data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point() 

# make it more tricky
data_loc = data_loc %>% slice(2:10)

# add grouping
point_n = nrow(data_loc)
group_n = round(point_n/2)
group_left = point_n%%2

point_identities = c( rep(0,each=group_left) , rep(1:group_n,each=2)  )
point_identities = point_identities + 1
data_loc = data_loc %>% mutate( point_group = point_identities   )

data_loc %>% ggplot(aes(x=truth_date,y=value)) + geom_line() + geom_point(aes(col=as.factor(point_group))) 
data_loc %>% ggplot(aes(x=truth_date,y=value_log)) + geom_line() + geom_point(aes(col=as.factor(point_group))) 



# make stan list

# fit stan model

# plotting 