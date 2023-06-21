#############################################################################################################
#### Script for runninng ap-sbw model #######################################################################
#### 

rm(list=ls())
source("R/ap-sbw.r")

scn="scn0"; #scn="scn2", scn="scn3", scn="scn4", scn="scn5", scn="scn6", scn="scn7", scn="scn8"
nrun=2; 
time.step=1; #only 5 
time.horizon=80; #0:80 & n*5
custom.params=NA; 
rcp='rcp45'#'rcp85',NA
is.sbw=T; 
is.harvesting=T
is.harvloc=F
save.land=T; 
time.save=10;
out.path="outputs/test1_s0"


############################################ RUN ap-sbw() ##################################################

## 1 basic run
r = ap.sbw(scn="scn0", is.sbw=T, is.harvesting=T, is.harvloc=F, is.harvprem=F, custom.params=NA, rcp='rcp45', 
           nrun=1, time.step=1, time.horizon=80, save.land=F, time.save=5, out.path=NA )


## Save outputs
r = ap.sbw(scn="scn10", is.sbw=T, is.harvesting=T, is.harvloc=F, is.harvprem=F, custom.params=NA, rcp='rcp45', 
           nrun=1, time.step=1, time.horizon=80, save.land=T, time.save=1, out.path="outputs/provaharvest_scn10" )


## Change default parameters of the list
source("R/default.params.r")  
custom.params = default.params()
custom.params$year.ini = 1990
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, is.harvloc=F, is.harvprem=F, custom.params=custom.params, rcp='rcp45', 
           nrun=3, time.step=1, time.horizon=80, save.land=F, time.save=5, out.path=NA )
  
  
## Change values of an innput table, eg. soil.suitability of SAB
#*** incloure custom tables a la funció principal

data(default.tables)
soil.suitability = default.tables[["soil.suitability"]]
soil.suitability[soil.suitability$spp=="SAB",2:6] = 0
default.tables[["soil.suitability"]] = soil.suitability
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, is.harvloc=F, is.harvprem=F, custom.params=custom.params, rcp='rcp85', 
           nrun=3, time.step=5, time.horizon=80, save.land=F, time.save=5, out.path=NA )

  sbw.outbreak(custom.params = NULL, custom.tables = tbl, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 10, nrun = 1, save.land = FALSE, time.save=5, out.path = NA) 

