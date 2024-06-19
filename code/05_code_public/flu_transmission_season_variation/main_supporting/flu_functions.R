


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

rep_warning_wed = function(df_rep,ind_name){
  # reporting warning for the case of reporting on weekdays other than Wednesday
  if (nrow(df_rep)==0) return(invisible(NULL))
  shouldbe_wednesday=df_rep$date %>% weekdays() %>% table() %>% names()
  warn1 = (length(shouldbe_wednesday)>1); warn2 = shouldbe_wednesday[1]!="Wednesday"
  if (warn1|warn2) {pr=paste("Warning: some",ind_name,"reports on days other than Wed \n"); cat(red(pr))}
  return(invisible(NULL))
}

data_into_all_season = function(data,params,withforce=F){
  
  file_doesnot_exist = !file.exists("output/all_season.Rdata")
  if (file_doesnot_exist|withforce==T) {
    # initiate
    df_collect = list()
    df_i = 1
    
    country_short_input_v=data$epi$erviss_ili_ari %>% 
      filter(target==params$SIR_simple$target) %>% 
      pull(country_short) %>% unique() ; length(country_short_input_v)
    
    for (country_short_input_i in country_short_input_v) {
      start_year = data$epi$erviss_ili_ari %>% filter(country_short==country_short_input_i) %>% 
        pull(date) %>% min() %>% year() %>% as.numeric()
      while( start_year<=params$latest_start_year ) {
        
        season     = paste0(start_year,"/",start_year+1)
        start_date = ymd(paste0(start_year,params$season_start_monthday))
        end_date   = ymd(paste0(start_year+1,params$season_end_monthday))
        date_v = seq(from=start_date,to=end_date,by="day")
        date_v_wed = date_v[weekdays(date_v)=="Wednesday"]
        start_year = start_year+1 # do this up here due to next; statements
        
        ## ili/ari
        data$epi$erviss_ili_ari %>% 
          filter(country_short == country_short_input_i) %>% 
          select(-country_short) %>% 
          filter( date%in%date_v ) -> xinc_iliari; rep_warning_wed(xinc_iliari,"ili/ari")
        
        if ( nrow(xinc_iliari) == 0 ) next;
        # fill the date gaps
        crossing(target=c("ILIconsultationrate","ARIconsultationrate"), 
                 date=date_v_wed,
                 agegroup=c("age_00_04", "age_15_64", "age_05_14", "age_65_99", "age_total")
        ) %>% 
          left_join(  xinc_iliari,by = join_by(target,date,agegroup) ) %>% 
          fill(c("agegroup", "target"),.direction = "downup") -> xinc_iliari
        
        ## typing_sentinel
        data$epi$erviss_typing_sentinel %>% 
          filter(country_short == country_short_input_i,date%in%date_v) %>% 
          filter(pathogen=="Influenza",pathogensubtype=="total") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age,-yearweek) -> xtyping_sent
        # fill the date gaps
        if (nrow(xtyping_sent)>=1) crossing( date=date_v_wed,
                                             indicator=c("detections","positivity","tests"  ) ) %>% 
          left_join( xtyping_sent, by = join_by(date,indicator) ) %>% 
          fill( c("pathogentype", "pathogensubtype"),.direction = "downup" ) -> xtyping_sent; rep_warning_wed(xtyping_sent,"sent_typing")
        
        # typing_nonsentinel
        data$epi$erviss_typing_nonsentinel %>% 
          filter(country_short == country_short_input_i,date%in%date_v) %>%
          filter(pathogen=="Influenza",pathogensubtype=="total") %>% 
          select(-country_short,-survtype,-countryname,-pathogen,-age,-yearweek) -> xtyping_nonsent
        # fill the date gaps
        if (nrow(xtyping_nonsent)>=1) crossing( date=date_v_wed,
                                                indicator=c("detections","positivity","tests"  ) ) %>% 
          left_join( xtyping_nonsent, by = join_by(date,indicator) ) %>% 
          fill( c("pathogentype", "pathogensubtype"),.direction = "downup" ) -> xtyping_nonsent
        xtyping_nonsent %>% group_by(date) %>%  mutate(
          value=ifelse(indicator=="positivity",
                       value[indicator=="detections"]/value[indicator=="tests"],
                       value)
        ) %>% ungroup() -> xtyping_nonsent; rep_warning_wed(xtyping_nonsent,"nonsent_typing")
        
        ## combine sentinel and non-sentinel
        xtyping_sent    %>% rename(value_sent=value)   -> x1
        xtyping_nonsent %>% rename(value_nonsent=value)-> x2
        xtyping_combined=left_join(x1,x2,by = join_by(date, indicator, pathogentype,pathogensubtype)) %>%
          mutate(value_add_narm=replace_na(value_sent,0)+replace_na(value_nonsent,0)  ) %>% 
          group_by(date) %>%  mutate(
            value_add_narm=ifelse(indicator=="positivity",
                                  value_add_narm[indicator=="detections"]/value_add_narm[indicator=="tests"],
                                  value_add_narm)
          ) %>% ungroup()
        
        ## data quality measures
        ili_sum=xinc_iliari %>% filter(target=="ILIconsultationrate") %>% summarise(x=sum(value,na.rm=T)) %>% pull(x)
        ili_quality=xinc_iliari %>% filter(target=="ILIconsultationrate") %>% mutate(v_q= !is.na(value)&(value>0) ) %>% summarise(x=mean( v_q )) %>% pull(x)
        ari_sum=xinc_iliari %>% filter(target=="ARIconsultationrate") %>% summarise(x=sum(value,na.rm=T)) %>% pull(x)
        ari_quality=xinc_iliari %>% filter(target=="ARIconsultationrate") %>% mutate(v_q= !is.na(value)&(value>0) ) %>% summarise(x=mean( v_q )) %>% pull(x)
        
        ntests_sent = xtyping_sent %>% filter(indicator=="tests") %>% summarise(msum=sum(value,na.rm=T)) %>% pull(msum)
        tests_sentinel_quality = xtyping_sent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5,na.rm=T)) %>% pull(x)
        ntests_nonsent = xtyping_nonsent%>% filter(indicator=="tests")%>% summarise(msum=sum(value,na.rm=T)) %>% pull(msum)
        tests_nonsentinel_quality = xtyping_nonsent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5,na.rm=T)) %>% pull(x)
        
        ## plotting
        xinc_iliari %>% ggplot(aes(date,value))+geom_line()
        xtyping_sent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        xtyping_nonsent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        
        ## skipping 
        # if ( nrow(xinc_iliari) < 10 ) next;
        # if ( sum_inc < 300 ) next;
        
        ## printing
        pr=paste0("> Run: ",country_short_input_i," | ",season,
                 " | ili/ari: ",my_comma(ili_sum),"/",my_comma(ari_sum),
                 " | tests-sent/nonsent/combined:",
                 my_comma(ntests_sent),"/",my_comma(ntests_nonsent),"/",my_comma(ntests_sent+ntests_nonsent),
                 "\n"); cat(green(pr))
        
        # put data together
        df_collect[[df_i]] = tibble(
          country_short=country_short_input_i,
          season=season,
          #
          ili_sum=ili_sum,
          ili_quality=ili_quality,
          ari_sum=ari_sum,
          ari_quality=ari_quality,
          tests_sentinel=ntests_sent,
          tests_sentinel_quality=tests_sentinel_quality,
          tests_nonsentinel=ntests_nonsent,
          tests_nonsentinel_quality=tests_nonsentinel_quality,
          # nested dataframes
          nest(xinc_iliari) %>% rename(inc_iliari=data),
          nest(xtyping_sent) %>% rename(typing_sentinel=data),
          nest(xtyping_nonsent) %>% rename(typing_nonsentinel=data),
          nest(xtyping_combined) %>% rename(typing_combined=data)
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
