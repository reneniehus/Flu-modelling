# look at Nick's data for RespiCast

mod_log = read_csv(
  file=paste0(
    "./output/flu-forecast-hub/ECDC-norrsken_green/","2023-12-20","-ECDC-norrsken_green","ILI",".csv"),show_col_types=F )
mod_log %>% filter(location=="CZ",output_type_id==0.975,horizon==4)
mod_log %>% filter(location=="CZ") %>% pull(value) %>% max() # 729.6846

x = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/flu-forecast-hub/main/model-output/ECDC-norrsken_green/2023-12-20-ECDC-norrsken_green.csv")
x %>% filter(location=="CZ") %>% pull(value) %>% max() # 257.928

mod_log = read_csv(
  file=paste0(
    "./output/flu-forecast-hub/ECDC-norrsken_green/","2023-12-27","-ECDC-norrsken_green","ILI",".csv"),show_col_types=F )
mod_log %>% filter(location=="CZ") %>% pull(value) %>% max() # 154.884


# test data, thus values are all NA

#### we need ILI 
df1 %>% filter(indicator=="ILIconsultationrate") %>% group_by(countryname) %>% summarise(
  num_year_weeks = n_distinct(yearweek[age=="total"])
) %>% arrange(desc(countryname)) %>%  print(n=Inf) 
#### we need ARI 
df1 %>% filter(indicator=="ARIconsultationrate") %>% group_by(countryname) %>% summarise(
  num_year_weeks = n_distinct(yearweek[age=="total"])
) %>% arrange(desc(countryname)) %>%  print(n=Inf)



x = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/ensembles/data-processed/EuroCOVIDhub-ensemble_all/2023-12-18-EuroCOVIDhub-ensemble_all.csv")
x2 = x %>% select(-n_models)
g(x2)
write_csv(x2,file = "./output/ensemble_fix.csv")
(x2$forecast_date) %>% table()

x = read_csv("https://raw.githubusercontent.com/EU-ECDC/Respiratory_viruses_weekly_data/main/data/nonSentinelSeverity.csv")
g(x)
(x$age) %>% table()
(x$indicator) %>% table() %>% names() %>% dput()
library(Hmisc)
Hmisc::describe(x)
library(skimr)
skim(x)

######################
df1 = read_csv(file="./data/from Nick end Nov test respicast/ILIARIRates.csv")
g(df1)
df1$survtype %>% table() # only: primary care syndromic
df1$indicator %>% table() # ARIconsultationrate ILIconsultationrate
df1$age %>% table() #  0-4 15-64  5-14   65+ total 
viewsummary(df1)
Hmisc::describe(df1)


df2 = read_csv(file="./data/from Nick end Nov test respicast/nonSentinelSeverity.csv")
g(df2)
df2$survtype %>% table() # non-sentinel
df2$indicator %>% table() %>% names() %>% dput() # "deaths", "hospitaladmissions", "hospitalinpatients", "ICUadmissions", "ICUinpatients"
df2$age %>% table() # 4400  4879  4106  5177  7027  5101
Hmisc::describe(df2)


df3 = read_csv(file="./data/from Nick end Nov test respicast/nonSentinelTestsDetections.csv")
g(df3)
df3$survtype %>% table() # non-sentinel
df3$indicator %>% table() # detections      tests
df3$pathogen %>% table() # Influenza RSV SARS-CoV-2
df3$age %>% table() # 
Hmisc::describe(df3) # covid cases might be in here!

df3 %>% filter(yearweek=="2023-W46") %>% g()

df4 = read_csv(file="./data/from Nick end Nov test respicast/sentinelTestsDetectionsPositivity.csv")
g(df4)
df4$survtype %>% table() #primary care sentinel
df4$indicator %>% table() # detections positivity tests
df4$pathogen %>% table() # Influenza RSV SARS-CoV-2 
df4$age %>% table() # total
Hmisc::describe(df4)
df4 %>% filter(yearweek=="2023-W46") %>% group_by(countryname) %>% 
  summarise( n_total=n(),
             n_flu=sum(pathogen=="Influenza"),
             n_rsv=sum(pathogen=="RSV")) %>% arrange(desc(n_total)) %>% head(20)

df5 = read_csv(file="./data/from Nick end Nov test respicast/variants.csv")
g(df5)
df5$survtype %>% table() # SARS-CoV-2 variant
df5$datasource %>% table() # GISAID  TESSy

df5$age %>% table() # total
Hmisc::describe(df5)
