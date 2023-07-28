msm_degree_distribution = function( d_max, PLOT_FLAG=FALSE ){
  
  #From Table 5.8 in https://www.ecdc.europa.eu/sites/default/files/documents/European-MSM-internet-survey-2017-findings.pdf
  prob_d_i_reported <-c( NA, #1-10 
                11.9, #11-20 
                5.2, #21-30
                2.3, #31-40
                1.4 )/100 #41-50
  
  prob_d_i_reported[ 1 ] <- 1 - 23.3/100 - 4.3/100 - sum( prob_d_i_reported[2:length( prob_d_i_reported )] )
  prob_d_i_reported <- prob_d_i_reported/sum( prob_d_i_reported )
  
  coarse_FLAG <- TRUE
  
  if ( coarse_FLAG ){
    
    df_d_i <- tibble( log_d_i = log( 1:length( prob_d_i_reported ) ), 
                      log_prob_d_i = log( prob_d_i_reported ) )
    
    power_law_d_i <- lm( log_prob_d_i~log_d_i, df_d_i )  
    
    d_i <- 1:d_max
    
    log_prob_d_i_model <- log( d_i ) * power_law_d_i$coefficients[ 'log_d_i' ]
    
     
  }else{
    
    
    df_d_i <- tibble( log_d_i = log( 1:length( prob_d_i_reported ) ), 
                      log_prob_d_i = log( prob_d_i_reported ) )
    
    power_law_d_i <- lm( log_prob_d_i~log_d_i, df_d_i )  
    
    d_i <- 1:d_max
    
    log_prob_d_i_model <- predict( power_law_d_i, tibble( log_d_i=log( d_i ) ))
  }
  
  prob_d_i_model <- exp( log_prob_d_i_model )
  prob_d_i_model <- prob_d_i_model/sum( prob_d_i_model )
  
  mout <- tibble( d_i = d_i,
                  prob_d_i = prob_d_i_model )
  
  if ( PLOT_FLAG ){
    mout %<>% mutate( prob_d_i_reported = NA )
    mout %<>% mutate( power_law_d_i_reported = NA )
    
    mout$prob_d_i_reported[ between( mout$d_i, 1, length( prob_d_i_reported ) ) ] <- prob_d_i_reported
    mout$power_law_d_i_reported[ between( mout$d_i, 1, length( prob_d_i_reported ) ) ] <- exp( power_law_d_i$fitted.values )
    
    # plot of final degree distribution
    print(  mout %>% 
              ggplot( ) + 
              geom_line( aes(x=log10( d_i ), y=log10( prob_d_i_model ) ), size=1 ) + 
              geom_point(aes(x=log10( d_i ), y=log10( prob_d_i_reported )), size=2 ) +
              xlab( "log10( degree d )" ) +
              ylab( "log10( Pr[ degree d ] )") +
              theme_minimal()+
              theme(text = element_text(size = 20)) )
      
    mout %<>% select( -prob_d_i_reported , -power_law_d_i_reported ) 
  }
  
  
  return( mout )
}
