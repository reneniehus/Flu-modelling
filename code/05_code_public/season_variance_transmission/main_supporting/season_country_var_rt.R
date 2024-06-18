# run special analysis to example seasonal and country variability of flu transmission potential

# ---- |-Load the data ----
rt_df = read_csv(file="output/rt_season_country.csv")
# filter out pandemic seasons
rt_df = rt_df %>% filter(!season%in%c("2019/2020","2020/2021","2021/2022"))
df = rt_df

# ---- |-Quick look ----
# brief look at variability
((df$Rnull)/median(df$Rnull) )%>% quantile(probs = c(0.1,0.5,0.9)) # 0.9115603 1.0000000 1.1262625
# Mean per country
df %>%
  group_by(country_short) %>%
  summarise(mean = mean(Rnull))
# boxplot per season
boxplot(Rnull ~ season,
        data = df
) # indicated that it makes sense to remove 1 or 2 covid winters
# boxplot per country
boxplot(Rnull ~ country_short,
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
  Rnull = df$Rnull
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
med_and_quantiles(extract(fit, "Rnull_relative_country_sim")$Rnull_relative_country_sim)
# Seasonal relative effect
med_and_quantiles(extract(fit, "Rnull_relative_season_sim")$Rnull_relative_season_sim)

Rnull_var = gather_draws(fit,Rnull_relative_season_sim) 
write_csv(Rnull_var,file="output/Rnull_relative_season_sim.csv")

Rnull_var$.value %>% 
  quantile(prob=c(0.1,0.2,0.5,0.8,0.9))

precis(fit,pars="Rnull_relative_season_sim")

mcmc_areas(fit,pars=c("Rnull_relative_season_sim") ) + 
  geom_vline(xintercept = c(0.2,-0.15,-0.10,-0.05,0.05,0.10,0.15,0.20),linetype="dashed") +
  coord_cartesian(xlim = c(-0.3,+0.25))+
  scale_y_discrete(labels=c("Rnull_relative_season_sim"="Between season Rt")) 
 
