SIR = function( A=NULL, 
                N_i=NULL, 
                delta=NULL, 
                beta=NULL, 
                R_0_target=NULL, 
                T_sample=1,  
                I_init=NULL,
                I_init_proportion=NULL, #scalar
                R_init_vector=NULL, #vector
                num_clusters = NULL,
                n=100, 
                fit_flag = FALSE ){
  
  if (sum( c( !is.null( delta ), !is.null( beta ), !is.null( R_0_target ) )) != 2 ){
    stop( "Specify exactly two of the parameters beta, delta, R_0_target")
  }
  
  if (sum( c( !is.null( I_init ), !is.null( I_init_proportion ) )) != 1 ){
    stop( "Specify exactly two of the parameters I_init, I_init_proportion")
  }
  
  if ( !is.null( num_clusters ) ){
    if ( dim( A )[1] > num_clusters ){
      #A is not necessarily symmetric!
      
      #spectral clustering, see, e.g., Algorithm 1 in https://projecteuclid.org/journals/annals-of-statistics/volume-48/issue-6/Clustering-in-Block-Markov-Chains/10.1214/19-AOS1939.short
      res_svd <- svds( A, num_clusters, nu = num_clusters, nv = num_clusters, opts=list( niter=1e4 ) )
      
      A_truncated <- res_svd$u %*% diag( res_svd$d ) %*% t( res_svd$v )
      
      res_kmeans <- kmeans( cbind( A_truncated, t( A_truncated ) ), num_clusters, iter.max = 1e2 )
      
      
      A_kmeans <- array( NA, dim=c( num_clusters, num_clusters ))
      R_init_kmeans <- array( NA, dim=c( num_clusters, 1 ))
      N_i_kmeans <- array( NA, dim=c( num_clusters, 1 ))
      for ( cluster_i in 1:num_clusters ){
        ind_i <- ( res_kmeans$cluster == cluster_i )
        
        N_i_kmeans[ cluster_i ] <- sum( N_i[ ind_i ] )
        
        if ( !is.null( R_init_vector )){
          R_init_kmeans[ cluster_i ] <- mean( R_init_vector[ ind_i ] )
        }
        
        for ( cluster_j in 1:num_clusters ){
          ind_j <- ( res_kmeans$cluster == cluster_j )
          
          A_kmeans[ cluster_i, cluster_j ] <- sum( A[ ind_i, ind_j ] )/sum( ind_i )
        }
      }
      
      N_i <- N_i_kmeans
      R_init_vector <- R_init_kmeans
      A <- A_kmeans
    }
  }
  
  
  #network size
  N <- dim( A )[1]
  
  if ( length( N_i ) != N ){
    stop( "Dimensions of N_i and A are inconsistent!" )
  }
  if (!is.null( I_init )){
    if ( length( I_init ) != N ){
      stop( "Dimensions of I_init and A are inconsistent!" )
    }
  }
  
  #the initial recovered state vector
  if ( is.null( R_init_vector )){
    R_init_vector <- array( 0, dim = c( N, 1 )) 
  }else{
    if ( any( dim( R_init_vector )  != c( N, 1 ) ) ){
      stop( "Dimensions of R_init_vector and A are inconsistent!" )
    }
  }
  
  
  #demography vector
  pi_vec <- N_i/sum( N_i )
  
  
  #obtain the principal eigenvector
  W_0 <- A #proportional to next generation matrix  ASSUMING HOMOGENEOUS spreading rates beta, delta
  res <- eigs( W_0, k=1, which="LM")
  x1 <- res$vectors #principal eigenvector
  
  if ( max( Im( x1 ) ) > 0 ){
    stop( "Principal eigenvector x1 non-real, should not happen (Perron-Frobenius). Check irreducibility and non-negativity of next generation matrix W !")
  }
  
  x1 <- Re( x1 )
  
  #by convention: use non-negative x1 
  if ( x1[1]<0 ){
    x1 <- -x1
  }
  
  if ( length( unique( sign( x1 )))>1 ){
    browser()
    stop( "Principal eigenvector x1 does not have only positive entries, should not happen (Perron-Frobenius). Check irreducibility and non-negativity of next generation matrix W !")
  }
  
  # all-one vector
  u <- array( 1, dim = c( N, 1 )) 
  
  #preallocate
  S <- array( NA, dim = c( N, n ))
  I <- array( NA, dim = c( N, n ))
  R <- array( NA, dim = c( N, n ))
  
  #initial condition:
  
  if( is.null( I_init  )){
    #choosing the initial condition I[,1] proportional to the principal eigenvector x1 results in 
    #an (almost) strictly increasing incidence, at least for the SIS model. See Corollary 2 in 
    #https://www.nas.ewi.tudelft.nl/people/Piet/papers/IEEE_TNSE2019_viral_state_dynamics_discrete_time_NIMFA.pdf
    inner_product_I_init_and_x1 <- c( I_init_proportion/( t( x1 ) %*% pi_vec ) ) #then I_sum[ 1 ] == t( pi_vec ) %*% I[,1] == I_init_proportion
    I[,1] <- inner_product_I_init_and_x1 * x1
    
  }else{
    I[,1] <- I_init
  }
  R[,1] <- c( R_init_vector )
  S[,1] <- u - I[,1] - R[,1]
  
  #if R_0_target is given as input, set either delta (if delta==NULL) 
  #or beta (if beta==NULL) such that R_0_target is attained
  if ( !is.null( R_0_target ) ){
    if ( is.null( delta )){
      #set curing rate delta such that R_0_target is attained
      delta <- 1
      W_0 <- beta/delta * diag( S[,1] ) %*% A #next generation matrix
      # W_0 <- beta/delta * A #next generation matrix
      res <- eigs( W_0, k=1, which="LM")
      R_0 <- Re( res$values[ 1 ] )
      delta <- delta * R_0/R_0_target
    } else if ( is.null( beta ) ){
      #set infection rate beta such that R_0_target is attained
      beta <- 1
      W_0 <- beta/delta * diag( S[,1] ) %*% A #next generation matrix
      # W_0 <- beta/delta * A #next generation matrix
      res <- eigs( W_0, k=1, which="LM")
      R_0 <- Re( res$values[ 1 ] )
      beta <- beta * R_0_target/R_0
    }
  }
  
  #infection rate matrix 
  B <- beta * A
  
  #discrete-time values
  B_T <- T_sample*B
  delta_T <- T_sample*delta
  
  if ( fit_flag ){
    
    R_0 <- NA
    
    I_sum <- rep( NA, n )
    
    I_old <- I[ , 1 ]
    R_old <- R[ , 1 ]
    S_old <- 1 - I_old - R_old 
    
    I_sum[ 1 ] <- sum( pi_vec * I_old )
    
    for ( k in 1:( n-1 ) ){
      I_new <- I_old - delta_T*I_old + S_old * c( B_T %*% I_old )
      R_new <- R_old + delta_T*I_old
      
      I_old <- I_new
      R_old <- R_new
      S_old <-  1 - I_old - R_old
      
      # I_sum[ k+1 ] <- t( pi_vec ) %*% I_new 
      I_sum[ k+1 ] <- sum( pi_vec * I_new ) 
    }
    
    R_sum <- sum( pi_vec * R[ , 1 ] ) + delta_T * cumsum( I_sum )
    S_sum <- 1 - I_sum - R_sum
    
  }else{
    
    
    #obtain R_0
    W_0 <- beta/delta * diag( S[,1] ) %*% A #next generation matrix
    res <- eigs( W_0, k=1, which="LM")
    R_0 <- Re( res$values[ 1 ] )
    
    for ( k in 1:( n-1 ) ){
      S_old <- S[ , k ]
      I_old <- I[ , k ]
      R_old <- R[ , k ]
      
      FOI <- diag( S_old ) %*% B_T %*% I_old
      
      S[ , k+1 ] <- S_old - FOI
      I[ , k+1 ] <- I_old - delta_T*I_old + FOI
      R[ , k+1 ] <- R_old + delta_T*I_old
    }
    
    S_sum <-  t( pi_vec ) %*% S 
    I_sum <-  t( pi_vec ) %*% I
    R_sum <-  t( pi_vec ) %*% R
    
  }
  
  mout <- list( t = ( 1:n ) * T_sample,
                S = S, 
                I = I, 
                R = R, 
                S_sum = c( S_sum ),
                I_sum = c( I_sum ),
                R_sum = c( R_sum ),
                beta = beta, 
                delta = delta, 
                R_0 = R_0,
                T_sample = T_sample,
                I_init_proportion = I_init_proportion,
                A = A )
  
  return( mout )
}
