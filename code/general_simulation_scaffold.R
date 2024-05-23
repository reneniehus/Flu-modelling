# ---- |-Simulations ----
#  inputs
#- for each axis a list with available IDs
#-

## select needed IDs within each list & create mapping of the IDs across axes -> axis_ids_simulate

## loop through axis_ids_simulate, 
for (sim_i in 1:nrow(axis_ids_simulate) ) {
  # prepare vaccination
  # prepare natural severity
  # prepare behaviour
  
  # run: data-to-transmission
  # run: transmission-to-severity 
  # run: combine all targets -> mysim
  
  # save
  axis_ids_simulate$sim[sim_i] = mysim
}
