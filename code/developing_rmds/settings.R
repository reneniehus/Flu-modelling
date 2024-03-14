# settings
p = list()
# parameters
p$ar = 60 # yearly attack rate
p$re = 0.8 #
p$ar_and_re = p$ar*p$re
p$Rt_for_ini_est = 4
p$R0 = 6

# discr_si(seq(from=0,to=10,length.out=40), mu = (3.0*1.2), sigma = (3.0*1.2^2)^(1/2) ) %>% plot()
p$si_mean=3.6 #  
p$si_sd=2.07

p$si_discr_dist = discr_si(seq(from=1,to=20,length.out=20), mu = p$si_mean, sigma = p$si_sd )


#
pop_myexample =120793614

# transformed parameters


# paths

# auxiliary data and others
aux = list()
# aux$pops = read_csv("https://raw.githubusercontent.com/european-modelling-hubs/covid19-forecast-hub-europe/main/data-locations/locations_eu.csv",
#                     show_col_types=F)
# write_csv(x = aux$pops, file="../Big data/eu_pops.csv")
aux$pops = read_csv("../../../Big data/eu_pops.csv")