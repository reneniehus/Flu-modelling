process_and_save = function(params=NULL, data=NULL, models_out=NULL,save_submission){
  pr=paste("Processing and saving ... \n"); cat(green(pr))
  
  ## ---- |-Put all countries together ----
  df_para = NULL
  df_submission = NULL
  df_data_summaries = NULL
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
    # observed burden
    x1 = models_out$mout[[i]]$season_ili_mean
    x_df = tibble(location=names(models_out$mout)[i],
                  season_ili_mean=x1)
    df_data_summaries = rbind(df_data_summaries,x_df)
  }
  
  ## ---- |-Sense checks ----
  # explore fits
  pdf(width = 12, height = 8, "code/03_report/fit_flip.pdf")
  for (i in 1:length(models_out$mout)) {
    # figures
    models_out$mout[[i]]$plot_fit_byage %>% print()
  }
  dev.off()
  
  
  # explore parameters
  ppar = list()
  df_para %>% filter(para=="prop_ili_mu") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) + 
    scale_x_log10() + labs(subtitle="prop_ili_mu") -> ppar$p1;ppar$p1
  df_para %>% filter(para=="SIR_ini_mu[1]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`))  + labs(subtitle="S_ini")-> ppar$p2; ppar$p2
  df_para %>% filter(para=="SIR_ini_mu[2]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10() + labs(subtitle="I_ini") -> ppar$p3; ppar$p3
  # explore submissions: ordering
  df_submission %>% group_by(scenario_id) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(cum_burden_log)
  df_submission %>% group_by(scenario_id,location) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(location,cum_burden_log) %>% View()
  #
  sum_ili = df_submission %>% filter(pop_group=="total_vaxTotal") %>% 
    group_by(scenario_id,location,horizon) %>% summarise(value=mean(value)) %>% ungroup() %>% 
    group_by(scenario_id,location) %>% 
    summarise(proj_burden=sum(value)) %>% 
    left_join( df_data_summaries, by="location" ) %>% 
    mutate( proj_burden_rel = ((proj_burden/season_ili_mean)*100) %>% round() ) %>% 
    arrange(location,scenario_id) 
  
  ## ---- |-Save ----
  if (save_submission) {
    library(arrow)
    sub_cols = c("round_id","scenario_id","target","location","pop_group",
                 "horizon","target_end_date","output_type","output_type_id","value")
    df_submission %>% 
      select(any_of(sub_cols) ) %>% filter(output_type_id%in%c(1:400)) %>% 
      mutate(horizon=as.integer(horizon),output_type_id=as.character(output_type_id)) %>% 
      filter(scenario_id%in%c("A","B","C","D","E","F")) %>% 
      filter(!is.na(horizon)) -> x
    el_distinct = c("round_id"=1,"scenario_id"=6,"target"=1,"locatoin"=18,
                    "pop_group"=15,"horizon"=43,"output_type_id"=400)
    cumprod(el_distinct)
    x %>% write_parquet("../Big data/2024_2025_1_FLU-ECDC-flumod.parquet")
  }
  
  
  # pack the list
  rep_list = list(
    N_countries_fit=models_out$mout %>% length(),
    sum_ili=sum_ili
    # plots
    
  )
  
  return(rep_list)
}