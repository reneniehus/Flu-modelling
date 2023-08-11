source("setup01.R")
source("DK_model/bayesian_functions.R")
# I_ini
wave_settings_SIR$I_ini = 0.0001 * 1/50;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$I_ini = 0.0001 * 1;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$I_ini = 0.0001
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean)) + labs(caption="I_ini")
# beta
wave_settings_SIR$beta = 0.2 * 0.8;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$beta = 0.2 * 1.2;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$beta = 0.2
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean)) + labs(caption="beta")
# S_ini
wave_settings_SIR$S_ini = 0.65;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$S_ini = 0.90;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$S_ini = 0.8
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean)) + labs(caption="S_ini")
# prop_severe
wave_settings_SIR$prop_severe = 0.02 * 10;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$prop_severe = 0.02 * 1;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$prop_severe = 0.02
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean)) + labs(caption="prop_severe")


# does the combination of I_ini and beta achieve the same as changing prop_severe?
# achieve different wave sizes with prop_severe
wave_settings_SIR$prop_severe = 0.02 * 5;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$prop_severe = 0.02 * 1;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$prop_severe = 0.02
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean))+ scale_y_log10() + labs(caption="prop_severe")
# achieve different wave sizes with beta
wave_settings_SIR$beta = 0.2 * 0.8;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$beta = 0.2 * 1.2;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$beta = 0.2
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean))+ scale_y_log10() + labs(caption="beta and I_ini")

# difference between beta effect and S_ini effect
# achieve different wave sizes with beta
wave_settings_SIR$beta = 0.2 * 0.8;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$beta = 0.2 * 1.2;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$beta = 0.2
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean))+ scale_y_log10() + labs(caption="beta and I_ini")
# achieving different wave sizes with S_ini
wave_settings_SIR$S_ini = 0.65;  wave_df1 = sim_wave_SIR(wave_settings_SIR,"low" )
wave_settings_SIR$S_ini = 0.90;    wave_df2 = sim_wave_SIR(wave_settings_SIR,"high")
wave_settings_SIR$S_ini = 0.8
bind_rows(wave_df1,wave_df2) %>% ggplot(aes(x=t_v,y=severe_obs,col=sim_name)) + geom_point(alpha=0.2) + 
  geom_point(alpha=0.2)+
  geom_line(aes(y=severe_mean))+ scale_y_log10() + labs(caption="S_ini")
