# case incidence - biggest caveat: strongly changing case ascertainment
path_data = "https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-truth/ECDC/truth_ECDC-Incident%20Cases.csv"
case_incidenc_EU = read_csv(path_data)
describe(case_incidenc_EU)
# filter for minimal working version
case_incidenc_EU = case_incidenc_EU %>% filter(location%in%c("DK",
                                                             "EE","HU","IT",
                                                             "LU","NL","NO","PT","SE"))
case_incidenc_EU = case_incidenc_EU %>% fill(value,.direction="down")
# for simple exmple: treat EU as a single epi region
case_incidenc_EU = case_incidenc_EU %>% group_by(date) %>% summarise(value=sum(value)) %>% 
  mutate(location="EU")

df_ii = case_incidenc_EU
