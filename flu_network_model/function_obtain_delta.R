obtain_delta = function( A = NULL,
                         tau = NULL,
                         N_i = NULL, 
                         T_sample = NULL,  
                         I_init = NULL,
                         I_init_proportion = NULL, #scalar
                         R_init_vector=NULL, #vector
                         num_clusters = NULL,
                         num_obs = NULL,
                         I_rep = NULL, 
                         delta_min = NULL, 
                         delta_max = NULL  ){
  
  delta_0 <- 0.9 * delta_min + 0.1 * delta_max
  
  res <- nloptr( x0 = delta_0,
                 eval_f = deviation_model,
                 lb = delta_min,
                 ub = delta_max,
                 opts = list("algorithm" = "NLOPT_LN_COBYLA",
                             "xtol_rel" = 1e-2,# 1e-8,
                             "maxeval" = 1e3 ),
                 A = A,
                 tau = tau,
                 N_i = N_i, 
                 T_sample = T_sample,  
                 I_init = I_init,
                 I_init_proportion = I_init_proportion, #scalar
                 R_init_vector = R_init_vector, #vector
                 num_clusters = num_clusters,
                 num_obs = num_obs,
                 I_rep = I_rep )
  
  delta_opt <- res$solution
  
  return( delta_opt )
}

deviation_model = function( x, 
                            A = NULL,
                            tau = NULL,
                            N_i = NULL, 
                            T_sample = NULL,  
                            I_init = NULL,
                            I_init_proportion = NULL, #scalar
                            R_init_vector=NULL, #vector
                            num_clusters = NULL,
                            num_obs = NULL,
                            I_rep = NULL){
  
  delta_this <- x
  beta_this <- delta_this * tau
  samples_per_day <- 1/T_sample
  
  res_SIR <- SIR( A = A,
                  N_i = N_i, 
                  delta = delta_this, 
                  beta = beta_this, 
                  T_sample = T_sample, 
                  I_init_proportion = I_init_proportion / delta_this ,
                  R_init_vector =  R_init_vector,
                  n = num_obs, 
                  fit_flag = TRUE  )
  
  #deviation to daily data
  deviation <- norm( res_SIR$R_sum[ seq( 1, num_obs, samples_per_day ) ] - res_SIR$R_sum[ 1 ] - cumsum( I_rep ), type="2" )
  
  return( deviation )
  
}