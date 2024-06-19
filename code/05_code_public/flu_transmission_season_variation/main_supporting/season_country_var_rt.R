# run special analysis to example seasonal and country variability of flu transmission potential

# ---- |-Load the data ----
rt_df = read_csv(file="output/rt_season_country.csv")
# filter out pandemic seasons
rt_df = rt_df %>% filter(!season%in%c("2019/2020","2020/2021","2021/2022"))
df = rt_df

# ---- |-Quick look ----
# brief look at variability
((df$Rnull_eff)/median(df$Rnull_eff) )%>% quantile(probs = c(0.1,0.5,0.9)) # 0.9115603 1.0000000 1.1262625
# Mean per country
df %>%
  group_by(country_short) %>%
  summarise(mean = mean(Rnull_eff))
# boxplot per season
boxplot(Rnull_eff ~ season,
        data = df
) # indicated that it makes sense to remove 1 or 2 covid winters
# boxplot per country
boxplot(Rnull_eff ~ country_short,
        data = df
)

# ---- |-Prepare for stan and fit ----
# prep data list
stan_list <- list(
  N = nrow(df),
  N_seasons = n_distinct(df$season),
  N_countries = n_distinct(df$country_short),
  country = (df$country_short) %>% fct_inorder() %>% as.numeric(),
  season = (df$season) %>% fct_inorder() %>% as.numeric(),
  season_id = (df$season) %>% fct_inorder() %>% levels() %>% enframe(),
  N_baseline_seasons = 3,
  baseline_seasons = c(4,5,6) %>% as.integer(),
  Rnull_eff = df$Rnull_eff
)
# run stan model
file_fit = "output/season_country_var_rt_fit.Rdata"
if ( params$debug==T&file.exists(file_fit) ){
  load(file=file_fit)
} else {
  fit <- rstan::stan(
    file = "stan/season_country_var_rt.stan",
    chains = 4, thin = 4, iter = 2500,
    seed = 12, cores = getOption("mc.cores", 1L),
    control = list(
      # adapt_delta=0.9,
      # max_treedepth=14
    ),
    data = stan_list
  )
 
  
  save(fit,file=file_fit)
}


# ---- |-Look at model estimates ----
med_and_quantiles <- function(x) {
  tibble(
    q1 = quantile(x, .1),
    q2 = quantile(x, .2),
    q8 = quantile(x, .8),
    q9 = quantile(x, .9)
  )
}
# Country relative effect
med_and_quantiles(extract(fit, "Rnull_eff_relative_country_sim")$Rnull_relative_country_sim)
# Seasonal relative effect
med_and_quantiles(extract(fit, "Rnull_eff_relative_season_sim")$Rnull_relative_season_sim)

Rnull_eff_var = gather_draws(fit,Rnull_relative_season_sim) 
write_csv(Rnull_eff_var,file="output/Rnull_relative_season_sim.csv")

Rnull_eff_var$.value %>% 
  quantile(prob=c(0.1,0.2,0.5,0.8,0.9))

precis(fit,pars="Rnull_eff_relative_season_sim")

mcmc_areas(fit,pars=c("Rnull_eff_relative_season_sim","Rnull_eff_relative_seasonbaseline_sim") ) + 
  geom_vline(xintercept = c(0.2,-0.15,-0.10,-0.05,0.05,0.10,0.15,0.20),linetype="dashed") +
  coord_cartesian(xlim = c(-0.3,+0.25))
 
