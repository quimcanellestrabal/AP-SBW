############################################################################################################
#### Script for organizing and plotting ap-sbw model results ###############################################
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2)

#data
ap.sbw_results <- readRDS("outputs/test1/ap.sbw_results.rds")

SppByAgeClass = ap.sbw_results$SppByAgeClass
SuitabilityClass = ap.sbw_results$SuitabilityClass    
LandChange = ap.sbw_results$LandChange            
ArtificialPlanting = ap.sbw_results$ArtificialPlanting   
ArtificialPlanting.sm = ap.sbw_results$ArtificialPlanting.sm 
Cuts = ap.sbw_results$Cuts                 
Cuts.sm = ap.sbw_results$Cuts.sm               
SBWDefol = ap.sbw_results$SBWDefol              
SBWKill = ap.sbw_results$SBWKill 


scn.name = "test1"


############ Function to plot age class distribution per species over time ############
ageclass.spp = function(scn.name){
  
  ## Read abundance of species per age class and reclassify in 
  ## 'young', 'mature', and 'old'.
  ap.sbw_results <- readRDS(paste0("outputs/",scn.name,"/ap.sbw_results.rds"))

  spp.age <- ap.sbw_results$SppByAgeClass
  spp.age.year <- mutate(spp.age, ymo=ifelse(age.class %in% c("C10", "C30", "C50"), "young",
                                             ifelse(age.class %in% c("C70", "C90"), "mature", "old"))) %>% 
    group_by(run, year, spp, ymo) %>% summarise(area=sum(area)) %>% 
    group_by(year, spp, ymo) %>% summarise(area=mean(area)) 
  
  ## Plot and save
  ylab1 = expression(paste(10^3 %.% km^2))
  p = ggplot(spp.age.year, aes(x=year, y=area/10^3, group=ymo)) + geom_line(aes(colour=ymo)) +
    geom_point(aes(color=ymo)) + scale_color_brewer(palette="Dark2") + theme_classic() +
    facet_wrap(~ spp, ncol=5, nrow=2) +  scale_y_continuous(ylab1)
  ggsave(p, filename = paste0(getwd(), "/plots/ageclass.spp_", scn.name, ".png"),
         width = 10, height = 6)
}




############ Function to plot suitability per species and bioclim domain over time ############
suitclass.spp = function(scn.name){
  
  #Read SuitabilityClass
  ap.sbw_results <- readRDS(paste0("outputs/",scn.name,"/ap.sbw_results.rds"))
  suit.class <- ap.sbw_results$SuitabilityClass
  
  suit.class.year <- suit.class%>% group_by(year, potential.spp) %>% 
    summarise(poor=mean(poor), med=mean(med), good=mean(good)) %>% 
    gather("suit", "area", 3:5)
    
  suit.class.biodomain <- suit.class%>% group_by(year, bioclim.domain, potential.spp) %>% 
    summarise(poor=mean(poor), med=mean(med), good=mean(good))
  
  ## Plot and save
  p = ggplot(suit.class.year, aes(x=year, y=area, group=suit)) + geom_line(aes(colour=suit))+
    geom_point(aes(color=suit)) + scale_color_brewer(palette="Dark2") + theme_classic() +
    facet_wrap(~ potential.spp, ncol=4, nrow=2)
  ggsave(p, filename = paste0(getwd(), "/plots/suitclass.spp_", scn.name, ".png"),
         width = 10, height = 6)
  
}




############ Function to plot landscape changes per species over time ############
landch.spp = function(scn.name, is.map=F){
  
  #Read SuitabilityClass
  ap.sbw_results <- readRDS(paste0("outputs/",scn.name,"/ap.sbw_results.rds"))
  land.ch <- ap.sbw_results$LandChange
  
  land.ch.year <- drop_na(land.ch)%>% group_by(year, spp, trans) %>% 
    summarise(n_obs=length(trans))
  
  ## Plot and save
  p = ggplot(land.ch.year, aes(x=year, y=n_obs, group=trans)) + geom_line(aes(colour=trans))+
    geom_point(aes(color=trans)) + scale_color_brewer(palette="Dark2") + theme_classic() +
    facet_wrap(~ spp, ncol=5, nrow=2)
  ggsave(p, filename = paste0(getwd(), "/plots/landch.spp_", scn.name, ".png"),
         width = 10, height = 6)
  
  
  if(is.map){
    #Get coordinates
    load("data/landscape.rda")
    land.ch = drop_na(left_join(land.ch, landscape[,c(1:3)], by="cell.id"))
    
    ##Plot and save
    land.ch$trans2 = as.numeric(as.factor(land.ch$trans))
    png(paste0(getwd(), "/plots/landch.map_", scn.name, ".png"))
    plot(rasterFromXYZ(land.ch[,c(7,8,9)]), col=c("red", "blue", "green")) #!! Millorar el format del mapa
    dev.off()
  }
}



############ Function to plot change of species due to artificial planting over time ############
ap.spp = function(scn.name){
  
  ## Read abundance of species per age class and reclassify in 
  ## 'young', 'mature', and 'old'.
  ap.sbw_results <- readRDS(paste0("outputs/",scn.name,"/ap.sbw_results.rds"))
  
  ap <- ap.sbw_results$ArtificialPlanting
  ap.year <- ap %>% unite("trans", initial.spp:ap.spp, sep = "-", remove = F) %>%
    group_by(year, ap.spp, trans) %>% 
    summarise(n_obs=length(trans))
  
  ## Plot and save
  p = ggplot(ap.year, aes(x=year, y=n_obs, group=trans)) + geom_line(aes(colour=trans))+
    geom_point(aes(color=trans)) + theme_classic() +
    facet_wrap(~ ap.spp, ncol=3, nrow=2)
  ggsave(p, filename = paste0(getwd(), "/plots/apptrans.spp_", scn.name, ".png"),
         width = 10, height = 6)
}



############ Function to plot change of species due to artificial planting over time ############
harv.spp = function(scn.name){
  
  ## Read abundance of species per age class and reclassify in 
  ## 'young', 'mature', and 'old'.
  ap.sbw_results <- readRDS(paste0("outputs/",scn.name,"/ap.sbw_results.rds"))
  
  harv <- ap.sbw_results$Cuts
  harv.year <- harv %>% group_by(year, spp) %>% 
    summarise(n_obs=length(spp))
  
  ## Plot and save
  p = ggplot(harv.year, aes(x=year, y=n_obs, group=spp)) + geom_line(aes(colour=spp))+
    geom_point(aes(color=spp)) + scale_color_brewer(palette="Dark2") + theme_classic() +
  ggsave(p, filename = paste0(getwd(), "/plots/harv.spp_", scn.name, ".png"),
         width = 10, height = 6)
}







############ Function to plot maps ############
# scn.name = "test1"; r=1; t=5; variable="transition"
# variable="age.class"or "transition" or "new_spp"

mapping = function(scn.name, r=1, t=5, variable=NA){
  
  #Load
  landscape <- readRDS(paste0("outputs/",scn.name,"/landscape_",r,"run_",t,"t.rds"))
  load(file="data/mask.rda")
  
  ##Arrange table
  kk = landscape[,c(1,2,match(variable,names(landscape)))]
  kk[,c(3)] = as.numeric(as.factor(kk[,c(3)]))
  kk[is.na(kk)] <- 0
  
  
  #Plot and save
  map = mask
  map[!is.na(mask[])] = kk[,c(3)]
  # png(paste0(getwd(), "/plots/",variable,".map_", scn.name, ".png"))
  plot(map, col=heat.colors(5)[5:1], main=paste(variable, 2020+t))
  # dev.off()

}


