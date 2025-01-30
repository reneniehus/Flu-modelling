# ---- |-Clear ----
gc() # clear environment & memory

# ---- |-Set up ----
source("code/01_main_supporting/setup.R")

# ---- |-Look at logit distributions ----
prob_mean = 0.10
prob_upp = 0.25 

rbeta(3000,shape1=2,shape2=3) %>% dens()

rnorm( 3000,mean=logit(prob_mean), sd=abs(logit(prob_mean)-logit(prob_upp)) ) %>% dens()
rnorm( 3000,mean=logit(prob_mean), sd=abs(logit(prob_mean)-logit(prob_upp)) ) %>% inv_logit() %>% dens()
