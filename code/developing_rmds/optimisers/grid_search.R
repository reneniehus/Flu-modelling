# optimising immunity_boe function
# ---- |-grid search  ----

f = function(...) compute_Rt_without_immunity(...,
                                              df_Rt_value_pop,
                                              for_optim = T,
                                              p=p)
grid_i = 1
df_grid = list()
7^4
for (NE_and_dark_factor_i in seq(from=25,to=50,length.out=10) ) {
  for ( wane_imm_fract_weekly_1_log_i in seq(from=-8,to=-0.2,length.out=10) ) {
    for ( Imm_2_partial_i in seq(from=0.3,to=0.9,length.out=15) ) {
      for ( wane_imm_fract_weekly_2_log_i in seq(from=-18,to=-8,length.out=10) ) {
        
        #######################################
        target_i = f( par=c(NE_and_dark_factor_i,
                                                      wane_imm_fract_weekly_1_log_i,
                                                      Imm_2_partial_i,
                                                      wane_imm_fract_weekly_2_log_i) 
                                                )
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
df_grid_all %>% ggplot(aes(x=NE_and_dark_factor,y=target)) + geom_point() + scale_y_log10()
df_grid_all %>% ggplot(aes(x=wane_imm_fract_weekly_1_log,y=target)) + geom_point()+ scale_y_log10()
df_grid_all %>% ggplot(aes(x=Imm_2_partial,y=target)) + geom_point()+ scale_y_log10()
df_grid_all %>% ggplot(aes(x=wane_imm_fract_weekly_2_log,y=target)) + geom_point()+ scale_y_log10()

df_min = df_grid_all %>% filter(target==min(target)) 
df_Rt_without_imm = compute_Rt_without_immunity(par=c(df_min$NE_and_dark_factor,
                                                      df_min$wane_imm_fract_weekly_1_log,
                                                      df_min$Imm_2_partial,
                                                      df_min$wane_imm_fract_weekly_2_log) , 
                                                df_Rt_value_pop , 
                                                for_optim=F,
                                                p=p_loop)

(par_min = c(df_min$NE_and_dark_factor,
  df_min$wane_imm_fract_weekly_1_log,
  df_min$Imm_2_partial,
  df_min$wane_imm_fract_weekly_2_log) %>% dput() %>% round(2) ) #  50.36  -0.40   0.67 -10.86

