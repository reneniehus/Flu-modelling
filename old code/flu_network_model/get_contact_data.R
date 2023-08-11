library("magrittr")
library( "socialmixr")

date_min <- as.Date("2022-01-01")

survey_list <- list_surveys()

survey_list %<>% filter( date_added >= date_min ) 


comix_v2_nl_survey$participants$country %>% unique()