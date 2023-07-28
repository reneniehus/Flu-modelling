library("magrittr")
library("dplyr")
library("tidyverse")
library("RSpectra")
library("zoo")
library("caTools")
library("nloptr")
library("profvis")
library("nleqslv")
library( "pracma")
library( "socialmixr")
library( "lsei")
library( "pheatmap" )
library( "fst" )
library( "readxl" )



path_core_functions <- "../../../../COVID-19-input/"

source( paste0( path_core_functions, "core_functions/function_obtain_demography.R" ))
source( "function_SIR_model.R" )
source( "function_final_size_vs_tau.R" )
source( "function_tau_given_final_size.R" )
source( "function_solution_complete_graph.R" )
source( "function_create_network.R" )
source( "function_obtain_delta.R" )
source( "msm_degree_distribution.R" )

select <- dplyr::select

set.seed( 1 )
contact_data <- "polymod"  
# contact_data <- "comix_2" 
country_i <- "Austria"
# country_i <- "Belgium"

t_start <- as.Date( "2022-03-01")

# Model parameters --------------------------------------------------------

#FIXME update values
N_pop <- 1e3 #number of all individuals
# N_pop <- 25.74e6 #number of all individuals
agegroups <- seq( 0, 70, 10 ) #minimum age of each agegroup
# pi_vec <- c( 0.2, 0.3, 0.5 ) #relative size of each group
# vacc_uptake <- c( 0, 0.4, 0.8 )
# VE <- 0.85


samples_per_day <- 30 #the larger, the slower the runtime. if chosen too small, then the range for the spreading rates may be too small (since, e.g., delta*T_sample < 1 to have a *stable* discrete-time model)
# d_max <- 50 #maximum degree
num_groups <- length( agegroups ) #number of groups

N_i <- obtain_demography( country_i, ten_year_brackets = TRUE )
N_i_AUS <- obtain_demography( "Australia", ten_year_brackets = TRUE )

#merge 70-79yr and 80+yr age group
N_i <- c( N_i[1:7], sum( N_i[ 8:9 ] ))
names( N_i )[8] <- "70+yr"

N_i_AUS <- c( N_i_AUS[1:7], sum( N_i_AUS[ 8:9 ] ))
names( N_i_AUS )[8] <- "70+yr"


N <- num_groups
dark_factor <- 3

# Load data ---------------------------------------------------------------
AUS_data <- read_csv( "./data/data_AUS_ILI.csv", show_col_types = FALSE )

AUS_data %<>% transmute( date = as.Date( date ), 
                         I_rep = value)

AUS_data %<>% filter( date >= t_start )

# AUS_data %>% ggplot( aes( x= date, y = I_rep )) + geom_line()

t_start <- AUS_data %>% pull( date ) %>% min( )
t_end <- AUS_data %>% pull( date ) %>% max( )

I_rep_sum <- sum( AUS_data$I_rep )

#interpolate NAs
AUS_data <- left_join( tibble( date=seq( t_start, t_end, by="day" )), AUS_data, by="date")
AUS_data %<>% arrange( date ) 
AUS_data %<>% mutate( I_rep = na.approx( I_rep ) )

#To ensure that the total number of cases is not affected by the interpolation above
AUS_data$I_rep <- I_rep_sum * AUS_data$I_rep /sum( AUS_data$I_rep )

# "Preprocess" data -------------------------------------------------------

#smooth data
# AUS_data %<>%  mutate( I_rep = runmean( I_rep, k = 21 ) )

#FIXME next two:
#from absolute to relative infections
AUS_data %<>%  mutate( I_rep = I_rep/N_pop )

#account for underascertainment
AUS_data %<>%  mutate( I_rep = dark_factor * I_rep )


#(optional) fill period after last reported case
post_length <- 1e3
if ( post_length>0 ){
  
  train_set_length <- 15*7
  a_lsei <- matrix( 1:train_set_length ) 
  a_lsei <- cbind( a_lsei, repmat( 1, nrow( a_lsei ), 1 ) )
  b_lsei <- tail( log( AUS_data$I_rep ), train_set_length )
  
  c_lsei <- a_lsei[nrow( a_lsei ),]
  d_lsei <- tail( b_lsei, 1 )
  
  e_lsei <- matrix( c( -1, 0 ), ncol=2 )
  f_lsei <- 0
  
  sol_lsei <- lsei(a = a_lsei, b = b_lsei, 
                   # c = c_lsei, d = d_lsei,
                   e = e_lsei, f = f_lsei)
  
  
  a_lsei_pred <- matrix( seq( train_set_length+1, train_set_length+post_length ) ) 
  a_lsei_pred <- cbind( a_lsei_pred, repmat( 1, nrow( a_lsei_pred ), 1 ) )
  val_pred <- exp( a_lsei_pred %*% sol_lsei )
  
  sum_I_future <- sum( val_pred )
  # plot( c( b_lsei, rep( NA, post_length)) )
  # points( x=1:train_set_length, y=a_lsei %*% sol_lsei, pch="x" )
  # points( x=seq( train_set_length+1, train_set_length+post_length ),
  #         y=val_pred, pch="." )
  
  # AUS_data %<>% bind_rows( tibble( date = seq( max( AUS_data$date )+ 1,
  #                                              max( AUS_data$date )+ post_length, by="days"),
  #                                  I_rep = val_pred ) )
  
  #######################
  #######################
  #######################
  # AUS_data %<>% bind_rows( tibble( date = seq( max( AUS_data$date )+ 1, 
  #                                               max( AUS_data$date )+ post_length, by="days"),
  #                                   I_rep = c(rep( NA, post_length-1 ), 0 ) ) ) 
  # 
  # AUS_data %<>% arrange( date )
  # 
  # AUS_data %<>% mutate( I_rep = na.approx( I_rep ))
  # AUS_data <- AUS_data[ 2:nrow(AUS_data), ]  
}

# Set up the network ------------------------------------------------------

A <- create_network( network_type = contact_data ,
                     country_i = "EU_average", 
                     agegroups = agegroups,
                     N_i = N_i_AUS )

# print( pheatmap( diag( N_i_AUS ) %*% A/max( diag( N_i_AUS ) %*% A ),
#                  display_numbers = T,
#                  main = paste0( "EU_average (normalised)" ),
#                  cluster_cols = FALSE,
#                  cluster_rows = FALSE ))



res <- eigs( A, k=1 )
x1 <- abs( Re( res$vectors ) )

#Set operatornorm to 1
A <- A/Re( res$values )


# Set up the model --------------------------------------------------------

#number of observations
n <- nrow( AUS_data ) * samples_per_day 

#sampling time
T_sample <- 1/samples_per_day 

# initial state
I_init <- x1/sum( x1 ) * AUS_data$I_rep[ 1 ]

#FIXME
R_init <- 5000 * I_init
# R_init <- vacc_uptake * VE

S_init <- 1 - I_init - R_init

# Plot final size versus tau ---------------------------------------------------

plot_final_size_vs_tau <- FALSE

if ( plot_final_size_vs_tau ){
  #obtain the final size as a function of the effective infection rate tau
  df_final_size <- final_size_vs_tau( A = A,
                                      S_init = S_init,
                                      R_init = R_init,
                                      N_i = N_i_AUS,
                                      tau_min = 0.5,
                                      tau_max = 2,
                                      num_tau = 1e2 )
  
  res_ev <- eigs( A, k=1 )
  print(
    df_final_size %>% 
      ggplot(aes( x=tau, y=final_size ))+
      geom_line(size=1)+ 
      xlab("Effective Infection Rate tau")+
      geom_vline(xintercept = 1/Re( res_ev$values ), linetype = "dashed", color="blue", size=1 )+
      annotate("text", x=1/Re( res_ev$values )*0.9, y=0.7, label="Ep. thr. 1/lambda_1", angle=90, size=7, color="blue")+
      # geom_vline(xintercept = 1, linetype = "dashed" )+
      # xlab("Basic Reproduction Number R_0")+
      ylim( 0,1 )+
      ylab("Final Size")+ 
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
      theme_bw()  +
      theme(text = element_text(size = 20)) 
  )
}


# Obtain tau from fitting the final size to data --------------------------

final_size <- sum( AUS_data$I_rep ) + sum_I_future + N_i_AUS %*% R_init /sum( N_i_AUS )

mout <- tau_given_final_size( A = A,
                              S_init = S_init,
                              R_init = R_init,
                              N_i = N_i_AUS,
                              final_size = final_size )

tau_opt <- mout$tau
R_0 <- mout$R_0
final_size_model <- mout$final_size

# Obtain delta from fitting to the incidence ------------------------------

delta_min <- 0.01
delta_max <- 1

delta_opt <- obtain_delta(  A = A, 
                            tau = tau_opt,
                            N_i = N_i_AUS,
                            T_sample = T_sample, 
                            I_init_proportion = AUS_data$I_rep[1],
                            R_init_vector = R_init, 
                            num_obs = n,
                            I_rep = AUS_data$I_rep, 
                            delta_min = delta_min, 
                            delta_max = delta_max )

beta_opt <- delta_opt * tau_opt


res_SIR_opt <- SIR( A = A, 
                    N_i = N_i_AUS, 
                    delta = delta_opt, 
                    beta = beta_opt, 
                    T_sample = T_sample,  
                    I_init_proportion = AUS_data$I_rep[1] / delta_opt,
                    R_init_vector = R_init, 
                    n = n, 
                    fit_flag = FALSE )


# Summarise model to daily data -------------------------------------------

R_model_daily <- res_SIR_opt$R_sum[seq( 1, n, samples_per_day )]  -  res_SIR_opt$R_sum[ 1 ]
I_model_daily <- diff( R_model_daily )
I_model_daily <- c( I_model_daily[1], I_model_daily )

# Plot fits ---------------------------------------------------------------


#Incidence I(t)
print( bind_rows( tibble( I = AUS_data$I_rep*N_pop, 
                          R = cumsum( AUS_data$I_rep*N_pop ), 
                          t = 1:length( AUS_data$I_rep ), 
                          type = "data" ),
                  tibble( I = I_model_daily*N_pop, 
                          R = R_model_daily*N_pop, 
                          t = 1:length( AUS_data$I_rep ), 
                          type = "model" ) ) %>% 
         ggplot( aes( x=t, y=I) ) +
         ggtitle( "Australia" ) +
         geom_line(aes(colour=type))+
         xlab( "Time t [in days]") +
         ylab( paste0( "Incidence I(t) per ", N_pop )) +
         theme_minimal()+
         theme(text = element_text(size = 20)) )

#Recovered R(t)
print( bind_rows( tibble( I = AUS_data$I_rep*N_pop,
                          R = cumsum( AUS_data$I_rep*N_pop ),
                          t = 1:length( AUS_data$I_rep ),
                          type = "data" ),
                  tibble( I = I_model_daily*N_pop,
                          R = R_model_daily*N_pop,
                          t = 1:length( AUS_data$I_rep ),
                          type = "model" ) ) %>%
         ggplot( aes( x=t, y=R) ) +
         ggtitle( "Australia" ) +
         geom_line(aes(colour=type))+
         xlab( "Time t [in days]") +
         ylab( "Recovered R(t) - R(0)") +
         theme_minimal()+
         theme(text = element_text(size = 20)) )
