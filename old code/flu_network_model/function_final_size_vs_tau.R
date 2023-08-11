#obtain the final size as a function of the effective infection rate tau
final_size_vs_tau = function( A = NULL,
                              S_init = NULL,
                              R_init = NULL,
                              N_i = NULL,
                              tau_min = 0.01,
                              tau_max = 5,
                              num_tau = 1e2){
  
  N <- dim( A )[1]
  pi_vec <- t( N_i/sum( N_i ))
  
  if (is.null( R_init )){
    R_init <- rep( 0, N ) 
  }
  
  tau_all <- seq( tau_min, tau_max, (tau_max-tau_min)/(num_tau-1) ) 
  R_inf_sum <- rep( NA, num_tau )
  
  
  mat_const <- rbind( diag( N ), -diag( N ) )
  vec_const <- rbind( array(0, dim = c( N,1 )), array(-1, dim = c( N,1 )))
  
  
  for ( counter in 1:length( tau_all ) ){
    tau_i <- tau_all[ counter ]
    
    theta0 <- matrix( rep( 0.5, N  ))
    
    res <- constrOptim( theta = theta0,
                        f = fn_deviation, 
                        grad = grad_fn_deviation,
                        ui = mat_const, 
                        ci = vec_const, 
                        mu = 1e-04, 
                        outer.iterations = 1e5, 
                        outer.eps = 1e-08,
                        S_init = matrix( S_init ),
                        R_init = matrix( R_init ),
                        tau = tau_i,
                        A = A )
    
    R_inf <- res$par
    
    R_inf_sum[counter] <- pi_vec %*% R_inf
  }
  
  
  res_lambda1 <- eigs( A, k=1 )
  lambda_1 <- Re( res_lambda1$values )
  
  return( tibble( final_size = R_inf_sum, 
                  tau = tau_all,
                  R_0 = tau_all * lambda_1) )
}


fn_deviation = function ( x, S_init, R_init, tau, A ){
  return( norm( 1 - x - S_init * exp( -tau * A %*% ( x - R_init ) ), type="2")**2 )
}


grad_fn_deviation = function ( x, S_init, R_init, tau, A ){
  xi <- S_init * exp( -tau*A %*% ( x - R_init ) )
  y <- 1 - x - xi

  df <- 2*tau*x* ( t( A ) %*% ( xi * y )) - 2 * y

  return( df )
}
