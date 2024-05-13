# run special analysis for rt seasonal and country variability

rt_df <- read_xlsx("code/special_analyses/rt_season_country.xlsx")

df <- rt_df %>% mutate(country = as.numeric(fct_inorder(country_short)),
                 season = as.numeric(fct_inorder(season)))




stan_list = list(
  N = nrow(df),
  N_counties = max(df$seasons),
  N_seasons = max(df$country),
  country = df$country,
  season = df$season,
  Rnull = df$Rnull
)
fit=rstan::stan(
  file='code/special_analyses/season_country_var_rt.stan',
  chains=8 ,thin=8,iter=300,
  seed=12, cores = getOption("mc.cores", 1L),
  control=list(
    #adapt_delta=0.9,
    #max_treedepth=14
  ),
  data=stan_list
) # X mins
#precis(fit02,pars = c("Rnull_eff"))







# Overall mean
rt_df <- rt_df %>% mutate(mean_val = mean(Rnull))

# Mean per country
rt_df %>% 
  group_by(country_short) %>% 
  summarise(mean = mean(Rnull))

# Mean per season
rt_df %>% 
  group_by(season) %>% 
  summarise(mean = mean(Rnull))

# Lm
# aov(Rnull ~ season, data = rt_df)

lm(Rnull ~ season  + country_short, data = rt_df, offset = mean_val)

boxplot(Rnull ~ season,
        data = rt_df)

boxplot(Rnull ~ country_short,
        data = rt_df)
