load_flu_data = function(params=NULL){
  # ---- |-Set up ----
  data = read_csv(file="./data/iliari_2015W40onwards.csv",show_col_types=F)
  # FIXME: load fresh data from the web and save locally
  
  data = data %>% 
    mutate(date=ISOweek2date(paste0(yearweek,"-3"))) %>% 
    mutate(age = case_when(
      age == "0-4" ~ "age_00_04",
      age == "5-14" ~ "age_05_14",
      age == "15-64" ~ "age_15_64",
      age == "65+" ~ "age_65_99",
      age == "total" ~ "age_total",
      age == "unk" ~ "age_unk",
    )) %>% 
    select( 
      country_short=countrycode,
      date=date, # see if we need to be more explicit
      agegroup_size=denominator,
      target=indicator,
      agegroup=age,
      value=value
      ) 
  
  return(data)
}


