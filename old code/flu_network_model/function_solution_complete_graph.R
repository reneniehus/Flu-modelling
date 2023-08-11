SIR_complete_graph_closed_form = function( S_init=NULL, 
                                           R_init= NULL, 
                                           I_init=NULL,
                                           beta=NULL, 
                                           N_i=NULL, 
                                           delta=NULL, 
                                           t=NULL ){
  
  # See section "Special Cases B" in:
  # Kermack, William Ogilvy, and Anderson G. McKendrick. 
  # "A contribution to the mathematical theory of epidemics." 
  # Proceedings of the Royal Society of London. Series A, 
  # Containing papers of a mathematical and physical character 115.772 (1927): 700-721.
  
  pi_vec <- N_i/sum( N_i )
  n <- length( t )
  N <- length( S_init )
  
  s_0 <- c( t( pi_vec ) %*% S_init )
  i_0 <- c( t( pi_vec ) %*% I_init )
  
  delta_T <- t[2] - t[ 1 ] # sampling time
  # beta <- beta * N
  
  sqrt_min_q <- sqrt( ( beta/delta*s_0 - 1 )**2 + 2*s_0*i_0*( beta/delta )**2 )
  Phi <- atanh( ( beta/delta*s_0 - 1 )/sqrt_min_q )
  
  gamma_t_apx <- ( delta/beta )**2/s_0*(  beta/delta*s_0 - 1 + sqrt_min_q * tanh( sqrt_min_q/2 * delta * t - Phi  )  )
  
  
  c_t_apx <- diff( gamma_t_apx )/delta_T/delta
  c_t_apx <- c( c_t_apx[1], c_t_apx )
  
  R_apx <- array( 1, dim=c( N, 1 )) %*% gamma_t_apx + matrix( R_init ) %*% rep( 1, n )
  I_apx <- array( 1, dim=c( N, 1 )) %*% c_t_apx
  
  R_sum <- t( pi_vec ) %*% R_apx 
  I_sum <- t( pi_vec ) %*% I_apx
  
  mout <- list( t = t,
                I = I_apx, 
                R = R_apx, 
                S_sum = 1 - c( I_sum ) - c( R_sum ),
                I_sum = c( I_sum ),
                R_sum = c( R_sum ) )
  
  return( mout )
}





