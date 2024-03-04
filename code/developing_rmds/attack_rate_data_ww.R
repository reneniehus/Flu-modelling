### more sources
# collection of WW data-observatory: https://wastewater-observatory.jrc.ec.europa.eu/media/bulletin_files/01_January_2024_Bulletin_250124.pdf

# use Danish wastewater 
# https://en.ssi.dk/covid-19/national-surveillance-of-sars-cov-2-in-wastewater
# there is a "lab" column
# eyeballing shows that this seems to have a muted impact
path_data="../Big input data/2024-01-17_dk_wastewater_data.csv"
x=read_csv(path_data)
g(x %>% filter(date>"2022-10-01",date<"2024-02-03") ); describe(x)
x %>% filter(date>"2022-10-01",date<"2024-02-03") %>% ggplot(aes(date,rna_mean_faeces)) + geom_line() # 
# ready for main
df_ii = x %>% select(date=date,
                     value=rna_mean_faeces) %>% 
  mutate(location="DK")


# compare 2 seasons
x %>% 
  summarise(
    burden_last_season = sum( rna_mean_faeces[date>"2022-09-01"&date<"2023-01-06"] ),
    burden_curr_season = sum( rna_mean_faeces[date>"2023-09-01"&date<"2024-01-06"] )
  ) %>% 
  mutate( curr_burden_rel = round(burden_curr_season/burden_last_season,2)  ) %>% 
  summarise(n_locations=n(),median_curr_burden_rel=median(curr_burden_rel)  ) # 2.73
###




### some feedback
# ROK: avoid first half of 2022 in this estimation basically


#### look into death burden last year versus this year
path_data = "https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-truth/ECDC/truncated_ECDC-Incident%20Deaths.csv"
x = read_csv(path_data)
clc()
g(x)
x %>% describe()
x %>% group_by(location) %>% 
  summarise(
    burden_last_season = sum( value[date>"2022-09-01"&date<"2023-01-27"] ),
    burden_curr_season = sum( value[date>"2023-09-01"&date<"2024-01-27"] )
  ) %>% 
  filter( !is.na(burden_last_season) & !is.na(burden_curr_season) ) %>% 
  mutate( curr_burden_rel = round(burden_curr_season/burden_last_season,2)  ) %>% 
  summarise(n_locations=n(),median_curr_burden_rel=median(curr_burden_rel)  )

##### look into the infection burden last year versus this year
# NL data downloaded here: https://www.infectieradar.nl/results
path_data="../Big input data/data_trendline.csv"
x=read.csv(path_data,sep = ";",dec=",")
g(x)
x %>% mutate(WEEK=ymd(WEEK)) %>% as_tibble() %>% filter(AGE_GROUP=="Totaal") -> x 
x %>% ggplot(aes(x=WEEK,y=ARI_CASE_PERC)) + geom_line()
x %>% 
  summarise(
    burden_last_season = sum( ARI_CASE_PERC[WEEK>"2022-09-01"&WEEK<"2023-01-06"] ),
    burden_curr_season = sum( ARI_CASE_PERC[WEEK>"2023-09-01"&WEEK<"2024-01-06"] )
  ) %>% 
  mutate( curr_burden_rel = round(burden_curr_season/burden_last_season,2)  ) %>% 
  summarise(n_locations=n(),median_curr_burden_rel=median(curr_burden_rel)  )



