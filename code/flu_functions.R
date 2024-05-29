


# Computing the severity factor due to vaccines, assuming no waning
vaccine_severity_nowane = function(
    vaccine_uptake, # [t,a] fraction of population vaccinated by time-bin and by age
    VE_severe # assumed reduction in severity by a typical administered dose 
){
  severity_factor_vaccines = 1 - ( vaccine_uptake*VE_severe )
  
  list_out = list(
    severity_factor_vaccines=severity_factor_vaccines
  )
  return(list_out)
}

# Computing the severity factor due to natural immunity
natural_severity = function(
    incident_infections, # [t,a] # infections by time-bin, and by age
    natual_severity_options=NULL # a list of options for natural severity module 
){
  severity_factor_natural = (incident_infections*0 + 1)
  
  list_out = list(
    severity_factor_natural=severity_factor_natural
  )
  return(list_out)
}

# Computing severe outcomes from infections
severity_factor = function(
    incident_infections, # [t,a] infections by time-bin, and by age
    severity_baseline, # [a] fraction of infections that is severe by age
    severity_factor_vaccines, # [t,a] modifying factor for severity due to vaccines
    severity_factor_natural, # [t,a] modifying factor for severity due to natural immunity
    severity_options=NULL # a list of options for severity module 
){
  # 1: combine the factors impacting severity
  severity_factor_combined = severity_factor_vaccines*severity_factor_natural
  # 2: modulate the baseline severity using the combined severity factors, get the effective severity
  severity_modulated = severity_factor_combined
  for (a_i in 1:length(severity_baseline) ) {
    severity_modulated[,a_i] = severity_baseline[a_i] * severity_factor_combined[,a_i]
  }
  # 3: compute the severe outcomes given infections and the effective severity 
  incident_severe = incident_infections*severity_modulated
  
  list_out = list(
    severity_factor_combined=severity_factor_combined, # [t,a] combined factors that impact the raw severity
    severity_modulated=severity_modulated, # [t,a] the effective severity after accounting for all severity factors
    incident_severe=incident_severe # [t,a] the severe indicator
  )
  return(list_out)
}

combine_all_targets_SIR_simple = function(date_v,
                                          incident_infections,
                                          vaccine_uptake,
                                          incident_severe) {
  mout=tibble(
    date=date_v,
    inc_infections=incident_infections[,1],
    inc_doses=vaccine_uptake[,1],
    inc_death=incident_severe[,1]
  )
  return(mout)
}

data_into_all_season = function(data,params,withforce=F){
  
  file_doesnot_exist = !file.exists("output/all_season.Rdata")
  if (file_doesnot_exist|withforce==T) {
    # initiate
    df_collect = list()
    df_i = 1
    
    country_short_input_v=data$erviss_ili_ari %>% 
      filter(target==params$SIR_simple$target) %>% 
      pull(country_short) %>% unique() ; length(country_short_input_v)
    
    for (country_short_input_i in country_short_input_v) {
      start_year = data$erviss_ili_ari %>% filter(country_short==country_short_input_i) %>% 
        pull(date) %>% min() %>% year() %>% as.numeric()
      while( start_year<=params$latest_start_year ) {
        
        season     = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,params$season_start_monthday))
        end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
        date_v_fit = seq(from=start_date,to=end_date,by="day")
        start_year = start_year+1 # do this up here due to next; statements
        
        ## ili
        data$erviss_ili_ari %>% 
          filter(country_short == country_short_input_i) %>% 
          select(-country_short) %>% 
          filter( date%in%date_v_fit ) -> xinc_iliari
        if ( nrow(xinc_iliari) == 0 ) next;
        # fill the date gaps
        crossing(target=c("ILIconsultationrate","ARIconsultationrate"), 
                 date=seq( min(xinc_iliari$date), max(xinc_iliari$date), by="week"),
                 agegroup=c("age_00_04", "age_15_64", "age_05_14", "age_65_99", "age_total")
        ) %>% 
          left_join(  xinc_iliari,by = join_by(target,date,agegroup) ) %>% 
          fill(c("agegroup", "target"),.direction = "downup") -> xinc_iliari
        
        ## typing_sentinel
        data$erviss_typing_sentinel %>% 
          filter(country_short == country_short_input_i,date%in%date_v_fit) %>% 
          filter(pathogen=="Influenza",pathogensubtype=="total") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age,-yearweek) -> xtyping_sent
        # fill the date gaps
        if (nrow(xtyping_sent)>=1) crossing( date=seq( min(xtyping_sent$date), max(xtyping_sent$date), by="week"),
                                             indicator=c("detections","positivity","tests"  ) ) %>% 
          left_join( xtyping_sent, by = join_by(date,indicator) ) %>% 
          fill( c("pathogentype", "pathogensubtype"),.direction = "downup" ) -> xtyping_sent
        ntests_sent = xtyping_sent %>% filter(indicator=="tests") %>% summarise(msum=sum(value)) %>% pull(msum)
        
        # typing_nonsentinel
        data$erviss_typing_nonsentinel %>% 
          filter(country_short == country_short_input_i,date%in%date_v_fit) %>%
          filter(pathogen=="Influenza",pathogensubtype=="total") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age,-yearweek) -> xtyping_nonsent
        # fill the date gaps
        if (nrow(xtyping_nonsent)>=1) crossing( date=seq( min(xtyping_nonsent$date), max(xtyping_nonsent$date), by="week"),
                                                indicator=c("detections","positivity","tests"  ) ) %>% 
          left_join( xtyping_nonsent, by = join_by(date,indicator) ) %>% 
          fill( c("pathogentype", "pathogensubtype"),.direction = "downup" ) -> xtyping_nonsent
        xtyping_nonsent %>% group_by(date) %>%  mutate(
          value=ifelse(indicator=="positivity",
                       value[indicator=="detections"]/value[indicator=="tests"],
                       value)
        ) %>% ungroup() -> xtyping_nonsent
        ntests_nonsent = xtyping_nonsent%>% filter(indicator=="tests")%>% summarise(msum=sum(value)) %>% pull(msum)
        
        ## data quality measures
        ili_quality = mean(xinc_iliari > 0) %>% round(2)
        tests_sentinel_quality = xtyping_sent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5)) %>% pull(x)
        tests_nonsentinel_quality = xtyping_nonsent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5)) %>% pull(x)
        sum_inc = sum(xinc_iliari$value) 
        
        
        ## plotting
        xinc_iliari %>% ggplot(aes(date,value))+geom_line()
        xtyping_sent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        xtyping_nonsent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        
        ## skipping 
        # if ( nrow(xinc_iliari) < 10 ) next;
        # if ( sum_inc < 300 ) next;
        
        ## printing
        pr=paste("> Running:",country_short_input_i,"| season:",season,
                 "| ili:",sum_inc,
                 "| test_sent:",ntests_sent,
                 "| test_nonsent:",ntests_nonsent,
                 "\n"); cat(green(pr))
        
        # put data together
        df_collect[[df_i]] = tibble(
          country_short=country_short_input_i,
          season=season,
          #
          ili_sum=sum_inc,
          ili_quality=ili_quality,
          tests_sentinel=ntests_sent,
          tests_sentinel_quality=tests_sentinel_quality,
          tests_nonsentinel=ntests_nonsent,
          tests_nonsentinel_quality=tests_nonsentinel_quality,
          # nested dataframes
          nest(xinc_iliari) %>% rename(inc_iliari=data),
          nest(xtyping_sent) %>% rename(typing_sentinel=data),
          nest(xtyping_nonsent) %>% rename(typing_nonsentinel=data)
        )
        df_i = df_i + 1
        
      } # season loop
    } # country loop
    all_season = bind_rows(df_collect)
    save(all_season,file="output/all_season.Rdata")
  }
  
  
  load(file="output/all_season.Rdata")
  return(all_season)
}
