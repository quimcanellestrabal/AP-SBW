######################################################################################
#### Script for preparing plots on qgis ##############################################
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(dplyr); library(rgdal)



#### Mapa amb les espècies ####

for (scn in c(10)){
  for (t in c(5)){
    for (y in seq(2025, 2100, 5)){
      land <- readRDS(paste0("outputs/test",t,"_scn",scn,"/landscape_1run_",y,"t.rds"))
      land$sppnum = as.numeric(land$spp)
      rr = rasterFromXYZ(land[,c(2,3,ncol(land))])
      
      if(!file.exists(paste0("proves/AnualSpp/test",t,"_scn",scn)))
        dir.create(file.path(paste0("proves/AnualSpp/test",t,"_scn",scn)), showWarnings = T) 
      
      writeRaster(rr, paste0("proves/AnualSpp/test",t,"_scn",scn,"/LandSpp_scn",scn,"_t",t,"_y",y,".tif"), format="GTiff")
      # 1=BOJ; 2=EPN; 3=ERS; 4=NonFor; 5=OTH.FEU.N; 6=OTH.FEU.S; 7=OTH.RES.N; 8=OTH.RES.S; 9=PET; 10=SAB 
      
      print( paste0("scn=",scn, "; test=",t, "; year=",y))
    }  
  }
}







#### Generar punts anuals ####

for (scn in c(10)){
  for (t in c(5)){
    for (y in seq(2025, 2100, 5)){
      landscape <- readRDS(paste0("outputs/test",t,"_scn",scn,"/landscape_1run_",y,"t.rds"))
      kk = filter(landscape[c(2,3,16)], landscape$ny.def>0)
      
      sp_points <- SpatialPoints(kk[,c(1,2)])
      sp_points$val <- kk$ny.def
      
      if(!file.exists(paste0("proves/AnualSpp/test",t,"_scn",scn)))
        dir.create(file.path(paste0("proves/AnualSpp/test",t,"_scn",scn)), showWarnings = T) 
      writeOGR(sp_points, dsn=paste0("proves/AnualSpp/test",t,"_scn",scn), layer=paste0("scn",scn,"_t",t,"_y",y), driver="ESRI Shapefile", overwrite_layer = T)
      
      print( paste0("scn=",scn, "; test=",t, "; year=",y))
      
    }
  }
}


