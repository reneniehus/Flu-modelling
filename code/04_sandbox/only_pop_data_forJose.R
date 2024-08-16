library(fst)
library(crayon)
library(tidyverse)
# initiate empty dataframe
pop_df = NULL
# 1: load the list of countries
xlocations = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/supporting-files/locations_iso2_codes.csv",show_col_types = F)
country_v = xlocations$location_name
# 2: loop through countries and add the pop data to dataframe
for (country_i in country_v) {
  pr=paste("> Loading pop data for:",country_i,"... \n"); cat(green(pr))
  read_file=paste0("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/miscellaneous/population/",country_i,"_aggr.csv")
  xdf = read_csv(read_file,show_col_types = FALSE)
  xdf$country = country_i
  pop_df=rbind(pop_df,xdf)
}
pop_df = pop_df %>% select(country,age_group,population)
