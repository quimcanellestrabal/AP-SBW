
#########################################################################################################
#### Script for organizing and plotting ap-sbw model maps ###############################################
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(dplyr); library(RColorBrewer); library(sf)


#testing colors
COL=heat.colors(10)[c(10,6,4,1)]
plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1))
rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)


# variable="ageClass"or variable="transition" or variable="new_spp" or variable="curr.intens.def"
scn.name = "test20_scn0"; r=1; t=2021; variable="curr.intens.def"

#Function
mapping = function(scn.name, r=1, t=2020, variable=NA, is.save=F, is.plot=F, is.shp=F){
  
  #Load
  landscape <- readRDS(paste0("outputs/",scn.name,"/landscape_",r,"run_",t,"t.rds"))
  load(file="data/mask.rda")
  
  ##Arrange table
  kk = landscape[,c(1,2,match(variable,names(landscape)))]
  nivells <- names(table(kk[,c(3)]))
  kk[,c(3)] = as.numeric(as.factor(kk[,c(3)]))
  kk[is.na(kk)] <- 0
  
  
  #Select color
  if(variable=="curr.intens.def"){COL=c("grey80",heat.colors(10)[c(6,4,1)])}
  if(variable=="ageClass"){
    COL=c("#FFFFBF",brewer.pal(6, "YlGn"))
    nivells=c("NA",nivells)}
  if(variable=="transition"){
    COL=c("#FFFFBF",brewer.pal(4, "Dark2"))
    nivells=c("NA",nivells)}
  if(variable=="new_spp"){
    COL=c("#FFFFBF",brewer.pal(10, "Paired"))
    nivells=c("NA",nivells)}
  
  #Plot and save
  map = mask
  map[!is.na(mask[])] = kk[,c(3)]
  
  if(is.save){
    
    png(paste0(getwd(), "/plots/maps/scn",iscn,"/",variable,"_", scn.name,"_run",r,"_t",t, ".png"))
    par(xpd = F)
    plot(map, col=COL, main=paste(variable, t, "r=",r), legend=F)
    par(xpd = TRUE)
    legend(x = 800000, y = 1000000, legend = nivells, fill = COL)
    dev.off()
  }
  
  if(is.plot){
    par(xpd = F)
    plot(map, col=COL, main=paste(variable, t, "r=",r), legend=F)
    par(xpd = TRUE)
    legend(x = 800000, y = 1000000, legend = nivells, fill = COL)
    
    if(is.shp){
      my_shp = st_read("C:/Users/Quim/OneDrive - Universitat de Girona/Dades/GIS/North America/PoliticalBoundaries_Shapefiles/NA_PoliticalDivisions/data/bound_l/boundary_l_v2.shp")
      plot(my_shp, add=T, color="black")
    }
  }
  
  
}
#exemple::mapping("test4",r=1,t=5, variable="transition",is.save=F, is.plot=T)


##LOOP mapping()
c.variable = c("ageClass","transition","new_spp","curr.intens.def") 

for (iscn in 0:8){
  cat(paste0("Plotting scn",iscn,"_test3"))
  for (ivar in 1:4){
    for (t in seq(2021, 2100, 1) ){
      print(t)
      mapping(scn.name=paste0("test20_scn",iscn),r=6,t=t, variable=c.variable[ivar],is.save=F, is.plot=T, is.shp=F)
    }
  }
}





