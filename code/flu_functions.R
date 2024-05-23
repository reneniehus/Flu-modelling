


# Computing the severity factor due to vaccines, assuming no waning
vaccine_severity_nowane = function(
    vaccine_uptake, # [t,a] fraction of population vaccinated by time-bin and by age
    VE_severe # assumed reduction in severity by a typical administered dose 
){
  severity_factor_vaccines = 1 - ( vaccine_uptake*VE_severe )
  
  list_out = list(
    severity_factor_vaccines
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
    severity_factor_natural
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
    severity_factor_combined, # [t,a] combined factors that impact the raw severity
    severity_modulated, # [t,a] the effective severity after accounting for all severity factors
    incident_severe # [t,a] the severe indicator
  )
  return(list_out)
}
