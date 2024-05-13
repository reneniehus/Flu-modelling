# run special analysis for rt seasonal and country variability

rt_df <- read_xlsx("code/special_analyses/rt_season_country.xlsx")

df <- rt_df %>% mutate(
  country = as.numeric(fct_inorder(country_short)),
  season = as.numeric(fct_inorder(season))
)

# prep data list
stan_list <- list(
  N = nrow(df),
  N_seasons = max(df$season),
  N_countries = max(df$country),
  country = df$country,
  season = df$season,
  Rnull = df$Rnull
)

# run stan model
fit <- rstan::stan(
  file = "code/special_analyses/season_country_var_rt.stan",
  chains = 8, thin = 8, iter = 300,
  seed = 12, cores = getOption("mc.cores", 1L),
  control = list(
    # adapt_delta=0.9,
    # max_treedepth=14
  ),
  data = stan_list
)

extract(fit, "a")$a %>% median()
extract(fit, "b")$b %>% median()

# Mean per country
rt_df %>%
  group_by(country_short) %>%
  summarise(mean = mean(Rnull))

# boxplot per season
boxplot(Rnull ~ season,
  data = rt_df
)

# boxplot per country
boxplot(Rnull ~ country_short,
  data = rt_df
)


med_and_quantiles <- function(x) {
  tibble(
    q1 = quantile(x, .1),
    q2 = quantile(x, .2),
    q8 = quantile(x, .8),
    q9 =  quantile(x, .9)
  )
}

# Country relative effect
med_and_quantiles(extract(fit, "Rnull_relative_country_sim")$Rnull_relative_country_sim)
# Seasonal relative effect
med_and_quantiles(extract(fit, "Rnull_relative_season_sim")$Rnull_relative_season_sim)
