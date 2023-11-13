# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Big section ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ---- |-Subsetion: More details ----

# libraries
start_time <- Sys.time()
library( caTools ) #for the function runmean
library( limSolve ) #for solving constrained linear least squares
library( sparsevar ) #for computing the spectral radius
library( glmc )
library( mgcv) # for splines
library( forecast ) # forecast
library(gamlss)
library( pracma )
library(magrittr) # better pipes
library(tidyverse) 
library(rstan)
library(scales)
library(ggpubr)
library(lubridate)
library(fitdistrplus)
library(readxl)
library(purrr)
library(tidybayes)
library(bayesplot)
library(rethinking)
library(patchwork)
library(viridis)
library(wrapr) # in wavefeature #
library(tidylog) # only temporarily
library(summarytools)
library(zoo)
library(here)
library(dagitty)
library(ISOweek)
library(EpiEstim)
library(fst)
library(tictoc)
library( "EcdcColors" )
library(cmdstanr)
library( scoringutils )
# libraries not to include 
# library(tsibble)

# remasking
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
extract <- rstan::extract
date <- lubridate::date
intersect <- base::intersect
setdiff <- base::setdiff
union <- base::union
rstudent <- rethinking::rstudent
expand <- tidyr::expand
map <- purrr::map
discard <- purrr::discard
col_factor <- readr::col_factor
combine <- gridExtra::combine
area <- patchwork::area
view <- summarytools::view
compare <- rethinking::compare
#
. %>% dfSummary %>% view() -> viewsummary
filter_log <- tidylog::filter
left_join_log <- tidylog::left_join
g = glimpse
detach(package:tidylog, unload = T)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Super basic functions ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clc = function() cat("\014")
mstand <- function(v){
  mout <- v/sum(v)
  return(mout)
}
# return the last element of a vector
mlast <- function(v ){
  mout = v[length(v)]
  return(mout)
}
# return the first element of a vector
mfirst <- function(v ){
  mout = v[1]
  return(mout)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Other Options ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rstan_options(auto_write = TRUE)
n_chains <-  2
options(mc.cores = parallel::detectCores())
# Reset R's most annoying default options
options(stringsAsFactors = FALSE, 
        scipen = 999, 
        dplyr.summarise.inform = FALSE,
        tibble.print_min=4)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Settings for plotting ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
theme_set(theme_gray() +
            theme(panel.grid.minor = element_blank()))
# amazing: overriding function defaults
ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Set3")
geom_interval <- function(...) ggdist::geom_interval(...,alpha=0.4)
geom_lineribbon <- function(...) ggdist::geom_lineribbon(...,alpha=0.4)
geom_ribbon <- function(...) ggplot2::geom_ribbon(...,alpha=0.4)

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Dark2")

mean_qi <- function(...) ggdist::mean_qi(...,.width=0.80)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Super basic helpers ##########
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# my own very basic functions
odds <- function(p){
  (p/(1-p)) -> mout
  return(mout)
}
inv_odds <- function( odds ){
  # body
  mp <- odds / (1+odds)
  return(mp)
} # try it: inv_odds( (0.4/0.6) )
odds_log <- function( p ) {
  log(p/(1-p)) -> mout
  return(mout)
}
logit <- odds_log
#
countries <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", 
               "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", 
               "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", 
               "Lithuania", "Luxembourg", "Malta", "Netherlands", "Norway", 
               "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", 
               "Sweden")
countries_short <- c("AT", "BE", "BG", "HR", "CY", "CZ", 
                     "DK", "EE", "FI", "FR", "DE", "GR", 
                     "HU", "IS", "IE", "IT", "LV", "LI", 
                     "LT", "LU", "MT", "NL", "NO", 
                     "PL", "PT", "RO", "SK", "SI", "ES", 
                     "SE")
# EL 
#EU_short("Greece") <- "EL"
# EU_short("Greece","EL")
EU_short <- function(name_long,greece="GR"){
  name_short = name_long
  for (i in 1:length(name_long)) name_short[i] <- countries_short[which(countries%in%name_long[i])]
  if (name_long=="Greece"&greece=="GR") name_short<-"GR"
  if (name_long=="Greece"&greece!="GR") name_short<-"EL"
  return(name_short)
}
#name_short=c("DE","PL","DE","PT")
EU_long <- function(name_short){
  name_long = name_short
  for (i in 1:length(name_long)) name_long[i] <- countries[which(countries_short%in%name_short[i])]
  return(name_long)
}

ecdc_weektodate <- function( year_week ){
  if ( any( nchar( year_week )!=7 ) ){
    stop( "Input must be of the format yyyy-ww !")
  }
  
  date_out <- ISOweek2date( paste0( substr( year_week, 1, 4 ), "-W", substr( year_week, 6, 7 ), "-1"  ) ) 
  
  return( date_out )
}

ecdc_datetoweek <- function( date_in ){
  if ( any( class( date_in )!="Date" ) ){
    stop( "Input must be a date !")
  }
  iso_week <- date2ISOweek( date_in )
  
  year_week <- paste0( substr( iso_week, 1, 4 ), "-", substr( iso_week, 7, 8 ) )
  return( year_week )
}

#  less simple functions, more specific to project
quantile_df <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

column_stats_ingroups = function( df , mycolumn,mygroup , ... ) {
  mycolumn = enquo(mycolumn)
  mygroup = enquo(mygroup)
  
  mysumm = df %>% ungroup() %>% 
    reframe( quantile_df( !!mycolumn , ... ), 
             .by = !!mygroup )
  return(mysumm)
}




# EU_long(c("DE","PL","DE","PT"))
# EU_short("Germany")
ggsave_as <- function(p,figname,height=10,width=16){
  ggsave(plot=p,filename=paste0(here(), "/figures/",figname,".pdf"),
         height=height,width=width,unit="cm")
  
}

ggsave_as_png <- function(p,figname,height=10,width=16){
  ggsave(plot=p,filename=paste0(here(), "/figures/",figname,".png"),
         height=height,width=width,unit="cm")
  
}

# very end: timing
end_time <- Sys.time()
cat("Run time setup01.R :", round(end_time - start_time,2) , "sec")

