## initiate model parameters ##
layer_mass <- rep(NA,t_extent); layer_mass[1] <- prod(0,temperature[1]) # g cm-2
layer_initial_mass <- rep(NA,t_extent); layer_initial_mass[1] <- layer_mass[1] # g cm-2
layer_remaining_mass <- rep(NA,t_extent); layer_remaining_mass[1] <- 1 #
layer_thickness <- rep(NA,t_extent); layer_thickness[1] <- layer_mass[1]/density # cm
layer_elevation <- rep(NA,t_extent); layer_elevation[1] <- layer_thickness[1] # cm
wt_height <- layer_elevation[1] # cm, required to calculate wt change
wet_proportion <- rep(NA,t_extent); wet_proportion[1] <- 1 #
layer_hydro_k <- rep(NA,t_extent); layer_hydro_k[1] <- k_param_a*(exp(k_param_b*layer_remaining_mass[1])) # cm yr-1
layer_transmissivity <- rep(NA,t_extent); layer_transmissivity[1] <- layer_thickness[1]*layer_hydro_k[1] # cm2 yr-1
peat_transmissivity <- layer_transmissivity[1] # cm2 yr-1, required to calculate wt change
time_counter <- 0 #
timestep <- 1/(annual_tsteps) # yr

## initiate model outputs ##
wtd_output <- transmit_output <- prod_output <- decay_output <- height_output <- mass_output <- rep(NA,t_extent)

## algorithm ##
for (i in 1:t_extent){
  print(paste("year",i))
  oxic_decay <- oxic_decay_base*Q10_oxic^((temperature[i]-base_temp)/10) # yr-1
  anoxic_decay <- anoxic_decay_base*Q10_anoxic^((temperature[i]-base_temp)/10) # yr-1
  wtd_sum <- 0 # cm, for calculating annual average
  repeat{ # repeat to calculate decay and water processes sub-annually
    time_counter <- time_counter+1
    if (time_counter > annual_tsteps) {break}
    for (j in 1:i){ # mass/thickness/elevation for each cohort
      layer_mass[j] <- layer_mass[j]*((1-wet_proportion[j])*exp(-timestep*layer_remaining_mass[j]^x_factor*oxic_decay)+
                                        wet_proportion[j]*exp(-timestep*layer_remaining_mass[j]^x_factor*anoxic_decay)) # g cm-2
      layer_remaining_mass[j] <- layer_mass[j]/layer_initial_mass[j] #
      layer_thickness[j] <- layer_mass[j]/density # cm
      if (j==1){
        layer_elevation[j] <- layer_thickness[j] # cm
      } else {
        layer_elevation[j] <- layer_thickness[j]+layer_elevation[j-1] # cm
      }
    }
    peat_height <- layer_elevation[i] # cm, peat height at top
    wt_height <- wt_height+timestep*dHdt(netprecip[i],porosity,peat_transmissivity,wt_height,lateral_extent) # cm, wt change
    if (wt_height <= 0) {print("wt_height warning"); wt_height <- 0.0000001} #
    if (wt_height > peat_height) {wt_height <- peat_height} # cm, standing water lost
    wtd_sum <- wtd_sum + (peat_height-wt_height) # cm, wtd summed up
    for (j in 1:i){ # wet proportion/conductivity/transmissivity for each cohort
      if (j==1){
        if (wt_height >= layer_elevation[j]) {wet_proportion[1] <- 1} else {wet_proportion[1] <- wt_height/layer_thickness[1]} #
      } else {
        if (wt_height >= layer_elevation[j]) {wet_proportion[j] <- 1} #
        else if (wt_height <= layer_elevation[j-1]) {wet_proportion[j] <- 0} #
        else {wet_proportion[j] <- (wt_height-layer_elevation[j-1])/layer_thickness[j]} #
      }
      layer_hydro_k[j] <- k_param_a*(exp(k_param_b*layer_remaining_mass[j])) # cm yr-1
      layer_transmissivity[j] <- layer_hydro_k[j]*layer_thickness[j]*wet_proportion[j] # cm2 yr-1
    }
    peat_transmissivity <- sum(layer_transmissivity[!is.na(layer_transmissivity)]) # cm2 yr-1
  }
  mass_after_decay <- sum(layer_mass[!is.na(layer_mass)]) # g cm-2, how much remained
  wtd_ave <- wtd_sum/annual_tsteps # cm, annual average wtd
  if (wtd_ave > 66.8){
    layer_mass[i+1] <- 0.0000001 #
  } else {
    layer_mass[i+1] <- prod(wtd_ave,temperature[i]) # g cm-2, new production annually
  }
  layer_initial_mass[i+1] <- layer_mass[i+1] # g cm-2
  layer_remaining_mass[i+1] <- 1 #
  wet_proportion[i+1] <- 0 #
  layer_thickness[i+1] <- layer_mass[i+1]/density # cm
  layer_elevation[i+1] <- peat_height+layer_thickness[i+1] # cm
  time_counter <- 0 #
  # output
  wtd_output[i] <- layer_elevation[i+1]-wt_height # cm
  transmit_output[i] <- peat_transmissivity # cm2 yr-1
  prod_output[i] <- layer_mass[i+1] # g cm-2
  if (i > 1){
    decay_output[i] <- mass_output[i-1]-mass_after_decay # g cm-2
  }
  height_output[i] <- layer_elevation[i+1] # cm
  mass_output[i] <- sum(layer_mass[!is.na(layer_mass)]) # g cm-2
}