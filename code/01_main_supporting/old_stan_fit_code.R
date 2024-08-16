if (F) {
  # the original stan sampling
  fit00=rstan::stan(
    file='./stan/SIR_multiseason_age_vax.stan',
    chains=1 ,thin=1,iter=300, # a "debug run"
    #chains=2, thin=2, iter=300, # a "long run" 
    seed=5, cores = getOption("mc.cores", 1L),
    control=list(
      # adapt_delta=0.95, # look into increasing this, 0.98 or 0.99
      max_treedepth=15 # look into increasing this to, 15, 20 ect
    ),
    data=stan_list #,init = init_fun
  ) # , run time will scale with iter, and is a function of adapt_delta and max_treedepth, and is a function of luck
  if (F) { # setting initial conditions for stan fit & checking fits
    # using a previous fit, set initial conditions
    ini_tune = F
    if (ini_tune==T) {
      load(path_fit) # loading fit00
      fit_means = get_posterior_mean(fit00) # extract mean estimates
      row.names(fit_means)[1:32] # print parameter names
      # mp="sigma_i";mcmc_areas(fit00,mp);precis(fit00,depth=3,mp)
      #
      myl = vector(mode = "list", length = 32)
      myl[1:32] = fit_means[1:32]
      names(myl) = row.names(fit_means)[1:32]
      init_fun = function(...) myl
      save(init_fun,file=paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata") )
      save_as_general = T
      if (save_as_general) save(init_fun,file="output/ini_SIR__multiseason_age_vax.Rdata")
    }
    if (ini_tune==F) {
      # a specific country run will look for a country-specfic initial contition function file, if that does not exist, it loads a general one
      if (file.exists(paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata"))){
        load(file=paste0("output/ini_SIR__multiseason_age_vax",target_input,country_short_input,".Rdata"))  
      } else{
        load(file="output/ini_SIR__multiseason_age_vax.Rdata")
      }
    }
    
  } else {
    load(path_fit)
  }
  # ---- |-Block to sense check etc, Extract parameters and plot, see convergence ----
  if (F){
    # plot without Rethinking package
    fit00@model_pars # see which parameters are there
    # see parameter names with dimensions
    fit_means = get_posterior_mean(fit00) # extract mean estimates
    row.names(fit_means)[1:17]
    fit_means[1:17,] %>% enframe()
    mp="SIR_ini[3,1,1]"; mcmc_areas(fit00,mp)
    mp="reciprocal_phi"; summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary
    
    mp="prop_ili_mu"; summary(fit00,pars=mp,probs = c(0.1, 0.9))$summary
    # https://rstudio.github.io/cheatsheets/bayesplot.pdf
    # for quick convergence check: n_eff Rhat ( see other packages than Rethinking)
    # using Rethinkingp akcage
    # precis(fit00,pars=c("SIR_ini_mu"),depth = 3)
    # precis(fit00,pars=c("Rnull_eff"),depth = 2)
    # precis(fit00,pars=c("prop_severe"),depth = 2)
    # precis(fit00,pars=c("sigma_i"),depth = 2) # 3.66
    # precis(fit00,pars=c("SIR_ini"),depth = 3) 
  }
  
}