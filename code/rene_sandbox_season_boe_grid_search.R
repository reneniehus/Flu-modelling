# ---- |-grid search  ----
grid_i = 1
df_grid = list()
6^3
for (NE_and_dark_factor_i in seq(from=1,to=40,length.out=1) ) {
  for ( wane_imm_fract_weekly_1_log_i in seq(from=-4.5,to=-0.001,length.out=10) ) {
    for ( Imm_2_partial_i in seq(from=0.10,to=0.6,length.out=12) ) {
      for ( wane_imm_fract_weekly_2_log_i in seq(from=-10,to=-4,length.out=7) ) {
        
        #######################################
        target_i = compute_Rt_without_immunity( par=c(wane_imm_fract_weekly_1_log_i,
                                                      Imm_2_partial_i,
                                                      wane_imm_fract_weekly_2_log_i) , 
                                                df_Rt_value_pop , 
                                                for_optim=T)
        df_grid[[grid_i]] = tibble(NE_and_dark_factor=NE_and_dark_factor_i,
                                   wane_imm_fract_weekly_1_log=wane_imm_fract_weekly_1_log_i,
                                   Imm_2_partial = Imm_2_partial_i,
                                   wane_imm_fract_weekly_2_log=wane_imm_fract_weekly_2_log_i,
                                   target=target_i)
        grid_i = grid_i + 1
        #########################################
      }
      
    }
  }
}
df_grid_all = bind_rows(df_grid)
# df_grid_all %>% ggplot(aes(x=NE_and_dark_factor,y=target)) + geom_point() + scale_y_log10()
df_grid_all %>% ggplot(aes(x=wane_imm_fract_weekly_1_log,y=target)) + geom_point()+ scale_y_log10()
df_grid_all %>% ggplot(aes(x=Imm_2_partial,y=target)) + geom_point()+ scale_y_log10()
df_grid_all %>% ggplot(aes(x=wane_imm_fract_weekly_2_log,y=target)) + geom_point()+ scale_y_log10()

df_min = df_grid_all %>% filter(target==min(target)) 
df_Rt_without_imm = compute_Rt_without_immunity(par=c(df_min$wane_imm_fract_weekly_1_log,
                                  df_min$Imm_2_partial,
                                  df_min$wane_imm_fract_weekly_2_log) , 
                            df_Rt_value_pop , 
                            for_optim=F)
# look 
df_Rt_without_imm %>% 
  filter(date>"2023-01-01") %>% 
  ggplot(aes(x=date,y=Imm)) + 
  geom_line(aes(y=Rt_without_immunity),col="red") + 
  geom_line() + 
  geom_hline(yintercept = 1) + facet_wrap(~location,scales="free_y") + 
  theme(axis.ticks.y=element_blank())
