##################### NÚRIA #####################
rm(list=ls())
setwd("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW")
devtools::load_all()


source("R/ap-sbw.r")

scn="scn0"; 
nrun=1; 
time.step=1; #only 5 
time.horizon=2; #0:80 & n*5
custom.params=NA; 
rcp='rcp45'#'rcp85',NA
is.sbw=T; 
is.harvesting=T
is.harvloc=F
save.land=T; 
time.save=10;
out.path="outputs/test1_s0"
params$current.duration = 0

