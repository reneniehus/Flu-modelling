create_network = function( network_type, 
                           N = NULL, 
                           network_parameters = NULL, 
                           country_i = NULL, 
                           agegroups = NULL, 
                           N_i = NULL ){
  
  if ( network_type=="ER_directed" ){
    
    A <- (array( runif( N*N ),dim = c( N, N )) <= network_parameters$p_ER  )
    
    A[ A==TRUE ] <- 1
    A[ A==FALSE ] <- 0
    
  }else if ( network_type=="ER_undirected" ){
    
    A <- array( runif( N*N ),dim = c( N, N ))
    A <- ( A + t( A ))/2
    
    A <- ( A <= network_parameters$p_ER  )
    
    A[ A==TRUE ] <- 1
    A[ A==FALSE ] <- 0
    
  }else if( network_type=="configuration_model"  ){
    #see https://math.stackexchange.com/questions/3676422/edge-probability-and-expected-number-of-edges-in-the-configuration-model
    d_i <- network_parameters$degree_sequence
    N_i <- network_parameters$size_clusters
    
    num_stubs <- 2*sum( N_i*d_i )
    prob_link_from_N_i_to_N_j <- ( d_i %*% t( d_i )) / ( num_stubs - 1 )
    
    #quotient matrix, same as B^pi in https://www.nas.ewi.tudelft.nl/people/Piet/papers/Chaos2021_Clustering_for_Epidemics.pdf (E.g., Theorem 1) 
    A <- prob_link_from_N_i_to_N_j %*% diag( N_i )
    
    #A_ij is the number of links from *all* nodes in cell N_j to *a single* node in cell N_i
    
  } else if ( network_type == "polymod" ){
    
    data( polymod )
    countries_polymod <- unique( polymod$participants$country )
    
    if ( country_i %in% countries_polymod ){
      res <- contact_matrix( polymod, 
                             countries = country_i, 
                             age.limits = agegroups, 
                             symmetric = TRUE)
      
      #the reduced-size contact matrix, i.e., A==\bar{B}/beta in eq. (36) of https://www.nas.ewi.tudelft.nl/people/Piet/papers/Chaos2021_Clustering_for_Epidemics.pdf
      A <- res$matrix
      
      #the proportion of individuals in each group
      N_i_times_A <- diag( res$demography$proportion ) %*% A
      
      A <- diag( 1/N_i ) %*% diag( res$demography$proportion ) %*% A
      
      #check if symmetric:
      if( norm(  N_i_times_A - t( N_i_times_A  ), type="2" ) / norm( diag( N_i_times_A ) , type="2" ) > 1e-10 ){
        stop( "Contact matrix not symmetric!")
      }
      
    }else{
      
      filename_matrix <- "./data/polymod_average.rds"
      
      if ( file.exists( filename_matrix )){
        
        N_i_times_A <- readRDS( file = filename_matrix )
        
        A <- diag( 1/N_i ) %*% N_i_times_A
        
      }else{
        N <- length( agegroups )
        N_i_times_A <- array( 0, dim=c( N, N ))
        
        for ( country_polymod in countries_polymod ){
          res <- contact_matrix( polymod, 
                                 countries = country_polymod, 
                                 age.limits = agegroups, 
                                 symmetric = TRUE)
          
          N_i_times_A <- N_i_times_A + diag(  res$demography$proportion ) %*% res$matrix
          
        }
        N_i_times_A <- N_i_times_A/length( countries_polymod )
        
        A <- diag( 1/N_i ) %*% N_i_times_A
        
        #check if symmetric:
        if( norm(  N_i_times_A - t( N_i_times_A  ), type="2" ) / norm( diag( N_i_times_A ) , type="2" ) > 1e-10 ){
          stop( "Contact matrix not symmetric!")
        }
        
        saveRDS( N_i_times_A, filename_matrix )
        
      }
      
    }
    
    
    
  }else if( network_type == "comix_2" ){
    
    
    countries_comix <- c( "Austria", 
                          "Belgium", 
                          "Denmark", 
                          "Croatia", 
                          "Estonia", 
                          "Greece", 
                          "Italy", 
                          "Poland", 
                          "Portugal")
    
    filename <- "./data/comix_v2.rds"
    
    if ( file.exists( filename )){
      survey_data <- readRDS( file = filename )
    }else{
      survey_data <- get_survey( "https://doi.org/10.5281/zenodo.7014556" )
      saveRDS( survey_data, filename )
    }
    
    if ( country_i %in% countries_comix ){
      res <- contact_matrix( survey_data, 
                             countries = country_i, 
                             age.limits = agegroups, 
                             symmetric = TRUE)
      
      #the reduced-size contact matrix, i.e., A==\bar{B}/beta in eq. (36) of https://www.nas.ewi.tudelft.nl/people/Piet/papers/Chaos2021_Clustering_for_Epidemics.pdf
      A <- res$matrix
      
      #the proportion of individuals in each group
      N_i_times_A <- diag( res$demography$proportion ) %*% A
      
      A <- diag( 1/N_i ) %*% diag( res$demography$proportion ) %*% A
      
    }else{
      
      N_i_times_A <- array( 0, dim=c( N, N ))
      
      for ( country_comix in countries_comix ){
        res <- contact_matrix( survey_data, 
                               countries = country_comix, 
                               age.limits = agegroups, 
                               symmetric = TRUE)
        
        N_i_times_A <- N_i_times_A + diag(  res$demography$proportion ) %*% res$matrix
        
      }
      N_i_times_A <- N_i_times_A/length( countries_comix )
      
      A <- diag( 1/N_i ) %*% N_i_times_A
    }
    
    
    #check if symmetric:
    if( norm(  N_i_times_A - t( N_i_times_A  ), type="2" ) / norm( diag( N_i_times_A ) , type="2" ) > 1e-10 ){
      stop( "Contact matrix not symmetric!")
    }
    
  }else{
    stop( "Unknown network type")
  }
  
  
  return( A )
}




