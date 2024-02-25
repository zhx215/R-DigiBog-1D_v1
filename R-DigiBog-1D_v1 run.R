### The R code for 1-D DigiBog model (in Fortran) presented in Morris et al. (2015) ###
### Author: Zhengyu Xia ###

## start ##
setwd("C:/Users/zhyxi/Dropbox/digibog")
rm(list=ls())
set.seed(888)

## parameter ##
t_extent <- 3000 # yr
annual_tsteps <- 365 #
lateral_extent <- 20000 # cm, radius

oxic_decay_base <- 0.042 # yr-1 (from Morris et al., 2015)
anoxic_decay_base <- 0.0002 # yr-1 (from Morris et al., 2015)
base_temp <- 6.29 # degree C (from Morris et al., 2015)
Q10_oxic <- 2.5 # (from Morris et al., 2015)
Q10_anoxic <- 2.5 # (from Morris et al., 2015)

density <- 0.1 # g cm-3 (from Morris et al., 2015)
porosity <- 0.3 # (from Morris et al., 2015)
k_param_a <- 31740 # cm yr-1, hydraulic conductivity model parameter (from Morris et al., 2015)
k_param_b <- 8 # hydraulic conductivity model parameter (from Morris et al., 2015)
drain_effi <- 2 # draining efficiency, circular = 2, ellipse = 1 
x_factor <- 0.5 # recalcitrance (from Morris et al., 2015)

## climate forcing data ##
temperature <- rnorm(t_extent,7,1) # degree C
netprecip <- rnorm(t_extent,100,15) # cm

## adding climate shift ##
precip <- netprecip
temp <- temperature

## function ##
prod <- function(z,t){ # water table depth cm, air temp degree C
  0.0001*(9.3+1.33*z-0.022*z^2)^2*(0.1575*t+0.009132) # g cm-2 yr-1
}
dHdt <- function(U,theta_d,K_sum,H,L){ # net precip cm/yr, porosity, hydraulic conductivity sum cm/yr, water table height cm, radius cm
  U/theta_d-(drain_effi*K_sum*H)/(L^2*theta_d)
}

## run model ##
source("R-DigiBog-1D_v1 algorithm.R")

## plot results ##
source("R-DigiBog-1D_v1 plot.R")
