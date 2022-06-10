#############################################################################################################
#### Script for runninng ap-sbw model #######################################################################
#### 



rm(list=ls())
source("R/ap-sbw.r")

scn="scn1"; #scn="scn2", scn="scn3", scn="scn4", scn="scn5", scn="scn6", scn="scn7", scn="scn8"
nrun=3; 
time.step=5; #only 5 
time.horizon=80; #0:80 & n*5
custom.params=NA; 
rcp='rcp85'#'rcp45',NA
is.sbw=F; 
is.harvesting=T
save.land=F; 
out.seq=NA; 
out.path=NA; 


############################################ RUN ap-sbw() ##################################################

## 1 basic run
r = ap.sbw(scn="scn3", is.sbw=F, is.harvesting=T, custom.params=NA, rcp='rcp85', 
           nrun=3, time.step=5, time.horizon=80, save.land=F, out.seq=NA, out.path=NA )


## Save outputs
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, custom.params=NA, rcp='rcp85', 
           nrun=2, time.step=5, time.horizon=40, save.land=T, out.seq=NA, out.path="outputs/test1" )


## Change default parameters of the list
source("R/default.params.r")  
custom.params = default.params()
custom.params$year.ini = 1990
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, custom.params=custom.params, rcp='rcp85', 
           nrun=3, time.step=5, time.horizon=80, save.land=F, out.seq=NA, out.path=NA )
  
  
## Change values of an innput table, eg. soil.suitability of SAB
#*** incloure custom tables a la funció principal

data(default.tables)
soil.suitability = default.tables[["soil.suitability"]]
soil.suitability[soil.suitability$spp=="SAB",2:6] = 0
default.tables[["soil.suitability"]] = soil.suitability
r = ap.sbw(scn="scn1", is.sbw=F, is.harvesting=T, custom.params=custom.params, rcp='rcp85', 
           nrun=3, time.step=5, time.horizon=80, save.land=F, out.seq=NA, out.path=NA )

  sbw.outbreak(custom.params = NULL, custom.tables = tbl, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 10, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

