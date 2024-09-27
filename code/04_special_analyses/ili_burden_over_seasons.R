# setting the seasons
season_other = c("2014/2015","2015/2016","2016/2017","2017/2018","2018/2019")
season_focus = c("2023/2024")

# get the data
all_season = data_into_all_season( data,params,withforce=F ); 
plot_all_season(all_season)

# compute the week-by-week ratio to the reference season
# then obtain the per country-season median of that ratio
# then obtain the per season EU/EEA median of that ratio
all_season %>% 
  select(country_short,season,respicompass_ili_plus) %>% 
  unnest(respicompass_ili_plus) %>% filter(agegroup=="age_65_99") %>% select(-agegroup) %>% 
  mutate(week=week(date)) %>% select(-date) %>% 
  filter(season=="2023/2024") %>% select(-season) %>% rename(value_ref=value) -> df_ref
all_season %>% 
  select(country_short,season,respicompass_ili_plus) %>% 
  unnest(respicompass_ili_plus) %>% filter(agegroup=="age_65_99") %>% select(-agegroup) %>% 
  mutate(week=week(date))  %>% select(-date) %>% 
  left_join_log(df_ref,by=c("country_short","week")) %>% mutate(val_rel = value/value_ref) %>% 
  #
  group_by(country_short,season) %>% summarise( val_rel=median(val_rel,na.rm=T) ) %>% 
  group_by(season) %>% summarise(ip_rel=median(val_rel,na.rm=T)) %>% 
  mutate(ip_rel=ifelse(!season%in%c(season_other,season_focus), 0,ip_rel ) )-> x
  

#
x  %>% ggplot(aes(x=season,y=ip_rel)) + 
  geom_bar(stat = "identity",fill="grey") + 
  geom_bar(data=. %>% filter(season%in%season_focus),stat = "identity",fill="black") +
  geom_hline(yintercept = 1,col="black",linetype="dashed") +
  labs(x="",y="EU/EEA median ILI plus burden realive to 2023/24") + theme_tidybayes() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p;p
ggsave_as_png(p,"flu_burden_regression",12*1.0,16*1.0)
