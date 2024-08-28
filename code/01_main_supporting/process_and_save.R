process_and_save = function(params=NULL, data=NULL, models_out=NULL){
  
  ## ---- |-Put all countries together ----
  df_para = NULL
  df_submission = NULL
  for (i in 1:length(models_out$mout)) {
    # parameter estimates
    x = models_out$mout[[i]]$pars_df
    x = as_tibble(x, rownames = "para")
    x$country = names(models_out$mout)[i]
    x$country_long = names(models_out$mout)[i] %>% EU_long()
    df_para = rbind(df_para,x)
    # submission df
    x = models_out$mout[[i]]$modelled_proj
    df_submission = rbind(df_submission,x)
    browser()
  }
  
  ## ---- |-Sense checks ----
  # explore parameters
  ppar = list()
  df_para %>% filter(para=="prop_ili_mu") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) + 
    scale_x_log10() -> ppar$p1
  df_para %>% filter(para=="SIR_ini_mu[1]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) -> ppar$p2
  df_para %>% filter(para=="SIR_ini_mu[2]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10()  -> ppar$p3
  # explore submissions: ordering
  df_submission %>% group_by(scenario_id) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(cum_burden_log)
  df_submission %>% group_by(scenario_id,location) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(location,cum_burden_log) %>% View()
  
  ## ---- |-Save ----
  library(arrow)
  sub_cols = c("round_id","scenario_id","target","location","pop_group",
               "horizon","target_end_date","output_type","output_type_id","value")
  df_submission %>% 
    select(any_of(sub_cols) ) %>% filter(output_type_id%in%c(1:200)) %>% 
    mutate(horizon=as.integer(horizon),output_type_id=as.character(output_type_id)) %>% 
    filter(scenario_id%in%c("A","B","C","D","E","F")) %>% 
    filter(!is.na(horizon)) %>% 
    write_parquet("../Big data/2024_2025_1_FLU-ECDC-flumod.parquet")
  
  # pack the list
  rep_list = list(
    N_countries_fit=models_out$mout %>% length(),
    # plots
    
  )
  
  return(rep_list)
}