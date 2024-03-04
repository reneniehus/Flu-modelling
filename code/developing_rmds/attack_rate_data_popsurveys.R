### some inspiration
# https://x.com/adamjkucharski/status/1746913088559808744?s=20

# https://www.rivm.nl/coronavirus-covid-19/onderzoek/infectieradar
# specifically here: https://dashboard.infectieradar.nl/data/data_positive_tests.csv
path_data="../Big input data/data_positive_tests.csv"
x=read.csv(path_data,sep = ";",dec = ",") %>% mutate(WEEK=ymd(WEEK))
x %>% select(date=WEEK,
             value=value) %>% 
  mutate(location="NL")
df_ii = x
###########################

g(x); describe(x)
x = x %>% 
  mutate(date=WEEK,value=POS_TEST_PERC)
x %>% filter(date>"2022-10-01",date<"2024-02-03") %>% ggplot(aes(date,value)) + geom_line() # 
# date, value, location

