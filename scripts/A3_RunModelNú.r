## Run the ap-sbw model in debug mode 
rm(list=ls())
setwd("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW")
devtools::load_all()

## Change some default parameters
params = default.params()
params$outbreak = 12   # so it will start with outbreak phase
params$w.host = 1

## Load initial conditions and clean the current outbreak from the landscape
load("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW/data/landscape.rda")
summary(landscape)
landscape$ny.def = 0
landscape$ny.def0 = 30
landscape$cum.intens.def = 0
landscape$curr.intens.def  = 0

## Run 
r = ap.sbw(scn="scn0", is.sbw=T, is.harvesting=T, is.harvloc=F, is.harvprem=F, custom.params=params, rcp='rcp45', 
           nrun=1, time.step=1, time.horizon=2, save.land=F, time.save=5, out.path="outputs/test", landscape = landscape)
  