rm(list=ls())

## Load the model in debug mode
setwd("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW")
devtools::load_all()

## Manually run the AP-SBW model --------------------------------------------------------------------------
## Change some default parameters
params = default.params()
params$preoutbreak = 2   # so it will start with the preoutbreak phase
# params$preoutbreak = 12   # so it will start with the outbreak phase
params$radius.outbreak.mid = 3
params$radius.outbreak.range = 1
params$time.horizon = 4
## Run it:
r = ap.sbw(scn="scn0", is.sbw=T, is.harvesting=T, custom.params=params, 
           rcp='rcp45', nrun=1, out.path="outputs/test") 
  

## Run the AP-SBW model for a set of scenarios --------------------------------------------------------------------------
## Read an excel file with a list of parameters values for each testing scenario,
## then, pass these values to the corresponding elements in the "params" list.
scenario.params = xlsx::read.xlsx("scripts/params_scenarios.xlsx", sheetName="scenario_params")
scenarios = c("scn01", "scn02"); scn = "scn01"
for(scn in scenarios){
  custom.params = default.params()  # to have a named list of the parameters 
  for(i in 1:length(custom.params))
    custom.params[[i]] = scenario.params[scenario.params$scenario==scn, names(custom.params)[i]]
  # transform class of the parameters
  custom.params$save.land.df = ifelse(custom.params$save.land.df %in% c("FALSE", "F"), F, T)
  custom.params$stop.end.phase = ifelse(custom.params$stop.end.phase %in% c("FALSE", "F"), F, T)
  custom.params$enable.succ = ifelse(custom.params$enable.succ %in% c("FALSE", "F"), F, T)
  custom.params$is.harvprem = ifelse(custom.params$is.harvprem %in% c("FALSE", "F"), F, T)
  # run the model
  out.path=paste0("outputs/", scn)
  res = ap.sbw(scn=scn, is.sbw=T, is.harvesting=T, custom.params=custom.params, rcp='rcp45', nrun=1, out.path)
  if(!file.exists(out.path))
    dir.create(file.path(out.path), showWarnings = T) 
  saveRDS(res, paste0(out.path, "/", "ap_sbw.rds"))
}


## Run the AP-SBW model from scratch --------------------------------------------------------------------------
## Load initial conditions and remove the current outbreak from the landscape to start from scratch
load("C:/Users/nuria.aquilue/OneDrive - ctfc.cat/QBCMOD/AP-SBW/AP-SBW/data/landscape.rda")
summary(landscape)
landscape$ny.def = 0
landscape$ny.def0 = 30
landscape$cum.intens.def = 0
landscape$curr.intens.def  = 0
## Run 
r = ap.sbw(scn="scn1", is.sbw=T, is.harvesting=T, is.harvprem=F, custom.params=params, rcp='rcp45', 
           nrun=2, out.path="outputs/test", landscape = landscape)