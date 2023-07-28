#obtain the effective infection rate tau from a given final size 
tau_given_final_size = function( A = NULL,
                                 S_init = NULL,
                                 R_init = NULL,
                                 N_i = NULL,
                                 final_size = NULL){
  
  N <- dim( A )[1]
  pi_vec <- t( N_i/sum( N_i ))
  
  if (is.null( R_init )){
    R_init <- rep( 0, N ) 
  }
  
  res_lambda1 <- eigs( A, k=1 )
  lambda_1 <- Re( res_lambda1$values )
  
  
  x0 <- matrix( c( rep( final_size/N, N  ), 0.5 ))
  x_lb <- array( 0, dim = c( N+1, 1 ))
  A_adj <- A
  A_eq <- c(  pi_vec, 0 )
  b_eq <- final_size
  
  res <-   fmincon( x0 = x0, 
                    fn = fn_deviation_R_inf_tau,
                    # gr = grad_fn_deviation_R_inf_tau, #...since gradient function is not used for SQP approach
                    S_init = matrix( S_init ),
                    R_init = matrix( R_init ),
                    A_adj = A_adj,
                    lb = x_lb,
                    Aeq = t( A_eq ), 
                    beq = b_eq, 
                    tol = 1e-10 )
  
  R_inf <- res$par[ 1:N ]
  tau <- res$par[ N+1 ]
  final_size_model <- c( pi_vec %*% R_inf ) 
  
  
  
  
  return( tibble( final_size = final_size_model, 
                  tau = tau,
                  R_0 = tau * lambda_1) )
}


fn_deviation_R_inf_tau = function ( x, S_init, R_init, A_adj ){
  
  N <- length( S_init )
  R_inf <- x[ 1:N ] 
  tau <- x[ N+1 ]
  
  return( norm( 1 - R_inf - S_init * exp( -tau * A_adj %*% ( R_inf - R_init ) ), type="2")**2 )
}


grad_fn_deviation_R_inf_tau = function ( x, S_init, R_init, A_adj ){
  N <- length( S_init )
  R_inf <- x[ 1:N ] 
  tau <- x[ N+1 ]
  
  xi <- S_init * exp( -tau*A_adj %*% ( R_inf - R_init ) )
  y <- 1 - R_inf - xi
  
  df_dR_inf <- 2*tau*R_inf* ( t( A_adj ) %*% ( xi * y )) - 2 * y
  
  df_dtau <- 2 * sum( y * xi * ( A_adj %*% R_inf ) )
  
  grad_f <- matrix( c( df_dR_inf, df_dtau ) )
  
  return( grad_f )
}
