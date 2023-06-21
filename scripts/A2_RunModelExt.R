#############################################################################################################
#### Script for runninng ap-sbw model in EXT COMPUTER########################################################
#### 

rm(list=ls())
source("R/ap-sbw.r")

test <- 20

for (iscn in c(0,1,2,3,4,5,6,7,8,9)){
  
  #Run model and save outputs
  ap.sbw(scn=paste0("scn",iscn), is.sbw=T, is.harvesting=T, is.harvprem=F, custom.params=NA, rcp='rcp45', 
     nrun=10, time.step=1, time.horizon=80, save.land=T, time.save=1, out.path=paste0("outputs/test",test,"_scn",iscn) )
  
}



#NoCC
for (iscn in c(0,1,2,3,4,5,6,7,8,9)){
  ap.sbw(scn=paste0("scn",iscn), is.sbw=T, is.harvesting=T, is.harvprem=F, custom.params=NA, rcp=NA, 
      nrun=10, time.step=1, time.horizon=80, save.land=T, time.save=1, out.path=paste0("outputs/test",test,"_scn",iscn,"_noCC") )
}

#PREM
for (iscn in c(0,1,2,3,4,5,6,7,8,9)){
  ap.sbw(scn=paste0("scn",iscn), is.sbw=T, is.harvesting=T, is.harvprem=F, custom.params=NA, rcp='rcp45', 
         nrun=10, time.step=1, time.horizon=80, save.land=T, time.save=1, out.path=paste0("outputs/test",test,"_scn",iscn,"_PREM") )
}
