# optimising immunity_boe function
# ---- |-bifurcation  ----

f = function(...) compute_Rt_without_immunity(...,
                                              df_Rt_value_pop,
                                              for_optim = T,
                                              p=p)
NE_and_dark_factor_max = 150
NE_and_dark_factor_min = 50

wane_imm_fract_weekly_1_max = -0.01
wane_imm_fract_weekly_1_min = -8

Imm_2_partial_max = 0.9
Imm_2_partial_min = 0.1

wane_imm_fract_weekly_2_max = -0.01
wane_imm_fract_weekly_2_min = -8

# Bisection in higher dimensions
# Assuming one optimum
# Doing bisection one dimension at a time

NE_and_dark_factor_mid_0 = mean(c(NE_and_dark_factor_min, NE_and_dark_factor_max))
wane_imm_fract_weekly_1_mid_0 = mean(c(wane_imm_fract_weekly_1_min, wane_imm_fract_weekly_1_max))
Imm_2_partial_mid_0 = mean(c(Imm_2_partial_min, Imm_2_partial_max))
wane_imm_fract_weekly_2_mid_0 = mean(c(wane_imm_fract_weekly_2_min, wane_imm_fract_weekly_2_max))

i=0
while (i < 10){
  print(i)
  ## NE_and_dark_factor ##
  # Mid points of the two intervals [min, mid0] and [mid0, max]
  NE_and_dark_factor_mid_a = mean(c(NE_and_dark_factor_min, NE_and_dark_factor_mid_0))
  NE_and_dark_factor_mid_b = mean(c(NE_and_dark_factor_mid_0, NE_and_dark_factor_max))
  eval_min = f(c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_a = f(c(NE_and_dark_factor_mid_a, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_b = f( c(NE_and_dark_factor_mid_b, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min) )
  eval_mid_0 = f( c(NE_and_dark_factor_mid_0, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_max = f( c(NE_and_dark_factor_max, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  # assuming you want to minimize abs(eval_x), where eval_x is a difference between two things (ie Rt[1] - Rt[2])
  if ((abs(eval_mid_a) < abs(eval_min)) & (abs(eval_mid_a) < abs(eval_mid_0))){
    NE_and_dark_factor_max = NE_and_dark_factor_mid_0
    NE_and_dark_factor_mid_0 = NE_and_dark_factor_mid_a
  } else if ((abs(eval_mid_b) < abs(eval_max)) & (abs(eval_mid_b) < abs(eval_mid_0))){
    NE_and_dark_factor_min = NE_and_dark_factor_mid_0
    NE_and_dark_factor_mid_0 = NE_and_dark_factor_mid_b
  } else if ((abs(eval_mid_a) > abs(eval_mid_0)) & (abs(eval_mid_b) > abs(eval_mid_0))){
    NE_and_dark_factor_min = NE_and_dark_factor_mid_a
    NE_and_dark_factor_max = NE_and_dark_factor_mid_b
  } else {
    stop("NE_and_dark_factor: either mid_0 value is the true minimum OR no minimum in range [min,max]")
  }
  
  ## wane_imm_fract_weekly_1 ##
  # Mid points of the two intervals [min, mid0] and [mid0, max]
  wane_imm_fract_weekly_1_mid_a = mean(c(wane_imm_fract_weekly_1_min, wane_imm_fract_weekly_1_mid_0))
  wane_imm_fract_weekly_1_mid_b = mean(c(wane_imm_fract_weekly_1_mid_0, wane_imm_fract_weekly_1_max))
  eval_min = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_a = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_mid_a, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_b = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_mid_b, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_0 = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_mid_0, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_max = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_max, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  # assuming you want to minimize abs(eval_x), where eval_x is a difference between two things (ie Rt[1] - Rt[2])
  if ((abs(eval_mid_a) < abs(eval_min)) & (abs(eval_mid_a) < abs(eval_mid_0))){
    wane_imm_fract_weekly_1_max = wane_imm_fract_weekly_1_mid_0
    wane_imm_fract_weekly_1_mid_0 = wane_imm_fract_weekly_1_mid_a
  } else if ((abs(eval_mid_b) < abs(eval_max)) & (abs(eval_mid_b) < abs(eval_mid_0))){
    wane_imm_fract_weekly_1_min = wane_imm_fract_weekly_1_mid_0
    wane_imm_fract_weekly_1_mid_0 = wane_imm_fract_weekly_1_mid_b
  } else if ((abs(eval_mid_a) > abs(eval_mid_0)) & (abs(eval_mid_b) > abs(eval_mid_0))){
    wane_imm_fract_weekly_1_min = wane_imm_fract_weekly_1_mid_a
    wane_imm_fract_weekly_1_max = wane_imm_fract_weekly_1_mid_b
  } else {
    stop("wane_imm_fract_weekly_1: either mid_0 value is the true minimum OR no minimum in range [min,max]")
  }
  ## Imm_2_partial ##
  # Mid points of the two intervals [min, mid0] and [mid0, max]
  Imm_2_partial_mid_a = mean(c(Imm_2_partial_min, Imm_2_partial_mid_0))
  Imm_2_partial_mid_b = mean(c(Imm_2_partial_mid_0, Imm_2_partial_max))
  eval_min = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_a = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_mid_a, wane_imm_fract_weekly_2_min))
  eval_mid_b = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_mid_b, wane_imm_fract_weekly_2_min))
  eval_mid_0 = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_mid_0, wane_imm_fract_weekly_2_min))
  eval_max = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_max, wane_imm_fract_weekly_2_min))
  # assuming you want to minimize abs(eval_x), where eval_x is a difference between two things (ie Rt[1] - Rt[2])
  if ((abs(eval_mid_a) < abs(eval_min)) & (abs(eval_mid_a) < abs(eval_mid_0))){
    Imm_2_partial_max = Imm_2_partial_mid_0
    Imm_2_partial_mid_0 = Imm_2_partial_mid_a
  } else if ((abs(eval_mid_b) < abs(eval_max)) & (abs(eval_mid_b) < abs(eval_mid_0))){
    Imm_2_partial_min = Imm_2_partial_mid_0
    Imm_2_partial_mid_0 = Imm_2_partial_mid_b
  } else if ((abs(eval_mid_a) > abs(eval_mid_0)) & (abs(eval_mid_b) > abs(eval_mid_0))){
    Imm_2_partial_min = Imm_2_partial_mid_a
    Imm_2_partial_max = Imm_2_partial_mid_b
  } else {
    stop("Imm_2_partial: either mid_0 value is the true minimum OR no minimum in range [min,max]")
  }
  
  ## wane_imm_fract_weekly_2 ##
  # Mid points of the two intervals [min, mid0] and [mid0, max]
  wane_imm_fract_weekly_2_mid_a = mean(c(wane_imm_fract_weekly_2_min, wane_imm_fract_weekly_2_mid_0))
  wane_imm_fract_weekly_2_mid_b = mean(c(wane_imm_fract_weekly_2_mid_0, wane_imm_fract_weekly_2_max))
  eval_min = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))
  eval_mid_a = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_mid_a))
  eval_mid_b = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_mid_b))
  eval_mid_0 = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_mid_0))
  eval_max = f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_max))
  # assuming you want to minimize abs(eval_x), where eval_x is a difference between two things (ie Rt[1] - Rt[2])
  if ((abs(eval_mid_a) < abs(eval_min)) & (abs(eval_mid_a) < abs(eval_mid_0))){
    wane_imm_fract_weekly_2_max = wane_imm_fract_weekly_2_mid_0
    wane_imm_fract_weekly_2_mid_0 = wane_imm_fract_weekly_2_mid_a
  } else if ((abs(eval_mid_b) < abs(eval_max)) & (abs(eval_mid_b) < abs(eval_mid_0))){
    wane_imm_fract_weekly_2_min = wane_imm_fract_weekly_2_mid_0
    wane_imm_fract_weekly_2_mid_0 = wane_imm_fract_weekly_2_mid_b
  } else if ((abs(eval_mid_a) > abs(eval_mid_0)) & (abs(eval_mid_b) > abs(eval_mid_0))){
    wane_imm_fract_weekly_2_min = wane_imm_fract_weekly_2_mid_a
    wane_imm_fract_weekly_2_max = wane_imm_fract_weekly_2_mid_b
  } else {
    stop("wane_imm_fract_weekly_2: either mid_0 value is the true minimum OR no minimum in range [min,max]")
  }
  
  if (abs(eval_min) < 0.01 ){
    stop("Done!")
  }
  i = i+1
}

f( c(NE_and_dark_factor_min, wane_imm_fract_weekly_1_min, Imm_2_partial_min, wane_imm_fract_weekly_2_min))