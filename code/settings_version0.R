settings = function() {
  params = list()
  # ---- |-Flu parameters ----
  params$models_to_run = c("SIR_simple") 
  params$rate_infectious = 0.2777778
  params$Rnull = 2.0; # https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-basic-reproductive-number-r0-for-epidemic-highly-pathogenic-avian-influenza-subtype-h5n1-spread/A60F72F5004F3BC5FAC2A3F8BB188A0F
  
  params$path_save_results = "output.fst"
  
  # ---- |-Data source ----
  
  # ---- |-Countries ----
  
  # ---- |-Fitting and uncertainty ----
  
  # ---- |-Simulation setting ---- 
  
  # ---- |-Flu scenerios ----
  return(params)
}

 



