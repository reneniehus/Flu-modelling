# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Big section ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## ---- |-Subsection: More details ----

start_time <- Sys.time()

## ---- |-Libraries ----
library(magrittr) # better pipes
library(tidyverse) # for all piping
library(readxl)
library(readr)
library(purrr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(viridis)
library(summarytools)
library(zoo)
library(here)
library(ISOweek)
library(scales)
library(fst)
library(crayon)
library(Hmisc)
library(lubridate)
library(data.table) # dataframes on steroids
library(dtplyr) # using data.table but with dplyr syntax

library(tidylog) # only temporarily

# remasking
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
date <- lubridate::date
intersect <- base::intersect
setdiff <- base::setdiff
union <- base::union
map <- purrr::map
discard <- purrr::discard
area <- patchwork::area
view <- summarytools::view
. %>% dfSummary %>% view() -> viewsummary
# keeping only select functions from tidy_log, then detach
filter_log <- tidylog::filter
left_join_log <- tidylog::left_join
fill_log <- tidylog::fill
detach(package:tidylog, unload = T)

# simplify calling functions
g = glimpse
