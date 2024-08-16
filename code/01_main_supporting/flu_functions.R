


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

mnaming = function(df,mnames){
  names(df) = mnames
  return(df)
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
    
    for (country_short_input_i in country_short_input_v) { # country_short_input_i = country_short_input_v[1] 
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
        
        ## respicompass ili_plus
        data$epi$respicompass_iliplus %>% 
          filter(country_short == country_short_input_i) %>% 
          select(-country_short) %>% 
          filter( date%in%date_v ) -> x_iliplus; rep_warning_wed(x_iliplus,"respicompass_iliplus")
        crossing(target=c("ili_plus"), 
                 date=date_v_wed,
                 agegroup=c("age_00_04", "age_15_64", "age_05_14", "age_65_99", "age_total")
        ) %>% 
          left_join(  x_iliplus,by = join_by(target,date,agegroup) ) %>% 
          fill(c("agegroup", "target"),.direction = "downup") -> x_iliplus
        
        
        ## data quality measures
        ili_sum=xinc_iliari %>% filter(target=="ILIconsultationrate") %>% summarise(x=sum(value,na.rm=T)) %>% pull(x)
        ili_quality=xinc_iliari %>% filter(target=="ILIconsultationrate") %>% mutate(v_q= !is.na(value)&(value>0) ) %>% summarise(x=mean( v_q )) %>% pull(x)
        ari_sum=xinc_iliari %>% filter(target=="ARIconsultationrate") %>% summarise(x=sum(value,na.rm=T)) %>% pull(x)
        ari_quality=xinc_iliari %>% filter(target=="ARIconsultationrate") %>% mutate(v_q= !is.na(value)&(value>0) ) %>% summarise(x=mean( v_q )) %>% pull(x)
        
        ntests_sent = xtyping_sent %>% filter(indicator=="tests") %>% summarise(msum=sum(value,na.rm=T)) %>% pull(msum)
        tests_sentinel_quality = xtyping_sent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5,na.rm=T)) %>% pull(x)
        ntests_nonsent = xtyping_nonsent%>% filter(indicator=="tests")%>% summarise(msum=sum(value,na.rm=T)) %>% pull(msum)
        tests_nonsentinel_quality = xtyping_nonsent %>% filter(indicator=="tests") %>% summarise(x=mean(value>5,na.rm=T)) %>% pull(x)
        
        ili_plus_sum=x_iliplus %>% summarise(x=sum(value,na.rm=T)) %>% pull(x)
        ili_plus_quality=x_iliplus %>% summarise(x=mean(!is.na(value))) %>% pull(x)
        
        ## plotting
        xinc_iliari %>% ggplot(aes(date,value))+geom_line()
        xtyping_sent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        xtyping_nonsent %>% filter(indicator=="positivity") %>% ggplot(aes(date,value))+geom_line()
        x_iliplus %>% filter(agegroup=="age_total") %>% ggplot(aes(date,value))+geom_line()
        
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
          ili_plus_sum=ili_plus_sum,
          ili_plus_quality=ili_plus_quality,
          # nested dataframes
          nest(xinc_iliari) %>% rename(inc_iliari=data),
          nest(xtyping_sent) %>% rename(typing_sentinel=data),
          nest(xtyping_nonsent) %>% rename(typing_nonsentinel=data),
          nest(xtyping_combined) %>% rename(typing_combined=data),
          nest(x_iliplus) %>% rename(respicompass_ili_plus=data)
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

transform_contracts = function(data,params) {
  
  #stop("Implement the 5th age group!")
  contacts_normalized_all = list()
  
  if (F){
    xlocations = read_csv(file="https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/supporting-files/locations_iso2_codes.csv",show_col_types = F)
    for (country_i in xlocations$location_name){ # country_i = xlocations$location_name[1]
      stop("Check if population matrices below neeed to be transposed or not")
      
      # Load original contact matrix
      contacts_orig = data$contact[[country_i]]
      if (length(contacts_orig) == 1 ){ # If contact matrix for this country is not available, skip and go to the next
        next
      }
      
      # Get population sizes
      read_file=paste0("https://raw.githubusercontent.com/european-modelling-hubs/RespiCompass/main/auxiliary-data/miscellaneous/population/",country_i,".csv")
      x_pop = read_csv(read_file,show_col_types = FALSE)
      x_pop_vec = x_pop$population 
      x_pop_vec = c(x_pop_vec[1:16], sum(x_pop_vec[17:21]))
      
      # Add the 80+ age group, assuming it has same per person nr of contacts as 75-79y age group
      x_contacts = cbind(contacts_orig, contacts_orig[,16])
      x_contacts = rbind(x_contacts, x_contacts[16,])
      contacts_orig = x_contacts
      
      # Fix the contact matrix non-symmetry issue by taking the mean value of the two (taking population size into account, obviously)
      contacts = NA*contacts_orig
      for (ii in 1:nrow(contacts_orig)){
        for (jj in 1:nrow(contacts_orig)){
          contacts[ii,jj] = (contacts_orig[ii,jj]*x_pop_vec[ii] + contacts_orig[jj,ii]*x_pop_vec[jj]) / (2*x_pop_vec[jj])
        }
      }

      # Get total number of contacts per age group; aka, each element is total number of contacts between age group i and j
      contacts_total = (as.matrix(contacts) * t(matrix(rep( x_pop_vec, 17), nrow = 17))) %>% round(digits = 1)
      # Need to use the transpose in the pop matrix above such that columns of the population matrix have the same element
      # This is because value contact[j,i] represents number of contacts of person in age group i with persons in age group j
      if (!isSymmetric(contacts_total, check.attributes = FALSE)){
        stop("The contacts_total matrix is not symmetric!")
      }
      
      ####
      #total_nr_contacts_per_person = sum( contacts_total[row(contacts_total)>=col(contacts_total)] ) / sum(x_pop_vec)
      #contacts = contacts_total / (t(matrix(rep( x_pop_vec, 17), nrow = 17)) * total_nr_contacts_per_person)
      ####
      
      # Change from 16 age groups to 4
      contacts_total_new = matrix(NA,4,4)
      #
      contacts_total_new[1,1] = contacts_total[1,1]
      #
      tmp_matrix = contacts_total[2:3,2:3] 
      contacts_total_new[2,2] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)]) # This part ensures that we add diagonal + only one off-diagonal part (e.g., upper but not lower), ensuring we don't count things twice
      contacts_total_new[1,2] = sum(contacts_total[1,2:3])
      contacts_total_new[2,1] = sum(contacts_total[2:3,1])
      #
      tmp_matrix = contacts_total[4:13,4:13]
      contacts_total_new[3,3] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)])
      contacts_total_new[1,3] = sum(contacts_total[1,4:13])
      contacts_total_new[3,1] = sum(contacts_total[4:13,1])
      contacts_total_new[2,3] = sum(contacts_total[2:3,4:13])
      contacts_total_new[3,2] = sum(contacts_total[4:13,2:3])
      #
      tmp_matrix = contacts_total[14:17,14:17]
      contacts_total_new[4,4] = sum(tmp_matrix[row(tmp_matrix)>=col(tmp_matrix)])
      contacts_total_new[1,4] = sum(contacts_total[1,14:17])
      contacts_total_new[4,1] = sum(contacts_total[14:17,1])
      contacts_total_new[2,4] = sum(contacts_total[2:3,14:17])
      contacts_total_new[4,2] = sum(contacts_total[14:17,2:3])
      contacts_total_new[3,4] = sum(contacts_total[4:13,14:17])
      contacts_total_new[4,3] = sum(contacts_total[14:17,4:13])
      #
      
      # Get total mean number of contacts per person
      total_nr_contacts_per_person = sum( contacts_total_new[row(contacts_total_new)>=col(contacts_total_new)] ) / sum(x_pop_vec)
      
      
      # Get a new contact matrix with only 4 age groups, such that average number of contacts per person equals to one
      x_new_pop = data$demography_respicast$population_pyramid %>% filter(country == country_i) %>% pull(population)
      x_pop_matrix = t(matrix(rep(x_new_pop,4), nrow=4))
      warning("Is the above above ok or should it be transposed?")
      
      # The new contact matrix where elements are per person contacts between age group i and j such that the (weighted) average number of contacts is 1
      contacts_normalized = contacts_total_new / (x_pop_matrix * total_nr_contacts_per_person)
      
      contacts_normalized_all[[country_i]] = contacts_normalized
      
      
      
    }
  }
  
  contacts_normalized_all = as_tibble( matrix(1/5,5,5) )
  
  return(contacts_normalized_all)
}

squash_axis <- function(from, to, factor) { 
  # A transformation function that squashes the range of [from, to] by factor on a given axis 
  
  # Args:
  #   from: left end of the axis
  #   to: right end of the axis
  #   factor: the compression factor of the range [from, to]
  #
  # Returns:
  #   A transformation called "squash_axis", which is capsulated by trans_new() function
  
  trans <- function(x) {    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(scales::trans_new("squash_axis", trans, inv, domain = c(from, to)))
}
