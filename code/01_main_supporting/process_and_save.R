process_and_save = function(params=NULL, models_out=NULL){
  
  ## ---- |-Put all countries together ----
  df_para = NULL
  df_submission = NULL
  for (i in 1:length(models_out$other)) {
    # parameter estimates
    x = models_out$other[[i]]$pars_df
    x = as_tibble(x, rownames = "para")
    x$country = names(models_out$other)[i]
    x$country_long = names(models_out$other)[i] %>% EU_long()
    df_para = rbind(df_para,x)
    # submission df
    x = models_out$other[[i]]$modelled_proj
    df_submission = rbind(df_submission,x)
  }
  
  ## ---- |-Sense checks ----
  # explore parameters
  df_para %>% filter(para=="prop_ili_mu") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) + 
    scale_x_log10()
  df_para %>% filter(para=="SIR_ini_mu[1]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`))
  df_para %>% filter(para=="SIR_ini_mu[2]") %>% 
    ggplot(  ) + 
    geom_pointrange(aes(y=country_long,x=mean,xmin=`10%`,xmax=`90%`)) +
    scale_x_log10()
  # explore submissions: ordering
  df_submission %>% group_by(scenario_id) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(cum_burden_log)
  df_submission %>% group_by(scenario_id,location) %>% 
    summarise(cum_burden_log=sum(value) %>% log()) %>% arrange(location,cum_burden_log) %>% View()
  
  ## ---- |-Save  ----
  library(arrow)
  df_submission %>% write_parquet("../Big data/2024_2025_1_FLU-ECDC-flumod.parquet")
  
  
  # Basic plot(s):
  # Show number of samples per model and date 
  mdf = NULL
  for (mod_id in params$models_to_run){
    mdf = bind_rows( mdf,models_out$df_for_submission[[mod_id]] )
  }
   
  avail_forecasts = 
    avail_forecasts(mdf %>% mutate(prediction=value, true_value=output_type_id, sample=output_type_id,model=model_id),
                    by=c("scenario_id","location","target_end_date","model"),
                    collapse = "sample_or_quantile")
  p = plot_avail_forecasts(
    avail_forecasts, x="target_end_date", show_numbers = F
  ) + facet_grid("scenario_tag"~"country_short")
  ggsave(filename = paste0(params$path_save_figures,"sample_overview.jpg"), plot = p)
  
}