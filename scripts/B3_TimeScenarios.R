#########################################################################################################-
# Script for organizing and analyzing ap-sbw model results ##############################################
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(dplyr)

#data
ap.sbw_results <- readRDS("outputs/test20_scn0/ap.sbw_results.rds")
landscape <- readRDS("outputs/test20_scn0/landscape_1run_2021t.rds")

SppByAgeClass = ap.sbw_results$SppByAgeClass
SuitabilityClass = ap.sbw_results$SuitabilityClass    
LandChange = ap.sbw_results$LandChange       
LandChange.sm = ap.sbw_results$LandChange.sm 
ArtificialPlanting = ap.sbw_results$ArtificialPlanting   
ArtificialPlanting.sm = ap.sbw_results$ArtificialPlanting.sm 
Cuts = ap.sbw_results$Cuts                 
Cuts.sm = ap.sbw_results$Cuts.sm               
SBWDefol = ap.sbw_results$SBWDefol 
SBWDefol.sm = ap.sbw_results$SBWDefol.sm   
SBWKill = ap.sbw_results$SBWKill
SBWKill.sm = ap.sbw_results$SBWKill.sm
SBWKill.sm2 = ap.sbw_results$SBWKill.sm2 
SppChange = ap.sbw_results$SppChange



### Cuts #####################################################################################################
Cuts = data.frame(run=NA, year=NA, a.age=NA, ncells=NA, scn=NA)
for (t in c(20)){
  for (iscn in c(1,2,3,4,5,6,7,8,9)){
    Cuts_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$Cuts.sm
    Cuts_scn$scn <- paste0("t",t,"_scn",iscn)
    Cuts <- rbind(Cuts, Cuts_scn)
  }
}

Cuts <- Cuts[-1,]

#Taula
taula <- Cuts %>%
  group_by(scn, year) %>%
  summarise(mean = mean(ncells), sd=sd(ncells))

#Plot
#taula = filter(taula, scn=="scn_8")
p1 = ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn)) + 
  geom_line() + geom_ribbon(alpha=0.5) +  xlab("Year") +  ylab("Harvested cells") 
p1
ggsave("plots/test20/Cuts.png", plot = p1, dpi = 300)



### LandChange ###############################################################################################
LandChange <- data.frame(run=NA, year=NA, trans=NA, n_obs=NA, scn=NA)

#load
for (t in c(20)){
  for (iscn in c(0:9)){
    LandChange_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$LandChange.sm
    LandChange_scn$scn <- paste0("t",t,"_scn",iscn)
    LandChange <- rbind(LandChange, LandChange_scn)
  }
}
LandChange <- LandChange[-1,]

#Taula
taula <- LandChange %>%
  group_by(scn, year, trans) %>%
  summarise(mean = mean(n_obs), sd=sd(n_obs))

#Plot
#taula = filter(taula, scn=="scn_8")
p2 = ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("Transition (ncell)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")+
  geom_ribbon(alpha=0.4, linetype = "blank") +
  scale_fill_brewer(palette = "Paired")+  
  facet_wrap(~ trans, ncol=5, nrow=1)
p2
ggsave("plots/test20/LandChange.png", plot = p2, dpi = 300)


p3 = ggplot(data=filter(taula, trans=="NatTrans"), aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("Transition (ncell)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")+
  geom_ribbon(alpha=0.4, linetype = "blank") +
  scale_fill_brewer(palette = "Paired")+  
  facet_wrap(~ trans, ncol=4, nrow=1)
p3
ggsave("plots/test20/LandChange_NatTrans.png", plot = p3, dpi = 300)


p4 = ggplot(data=filter(taula, trans=="Harv_AP"|trans=="Harv_ForTrans"), aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("Transition (ncell)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")+
  geom_ribbon(alpha=0.4, linetype = "blank") +
  scale_fill_brewer(palette = "Paired")+  
  facet_wrap(~ trans, ncol=4, nrow=1)
p4
ggsave("plots/test20/LandChange_Harv.png", plot = p4, dpi = 300)






### SBWKilled ###############################################################################################
SBWKill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA, scn=NA)

#load
for (t in c(20)){
  for (iscn in c(0:9)){
    SBWKill.sm_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$SBWKill.sm
    SBWKill.sm_scn$scn <- paste0("t",t,"_scn",iscn)
    SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
  }
}
SBWKill.sm <- SBWKill.sm[-1,]

#Taula
taula <- SBWKill.sm %>%
  group_by(scn, run, year) %>%
  summarise(km2 = sum(area)) %>%
  group_by(scn, year) %>%
  summarise(mean = mean(km2), sd=sd(km2))
taula[is.na(taula)] <- 0

#Plot
#taula = filter(taula, scn=="scn_8")
p5 = ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("SBWkilled (km2)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")#+
  #geom_ribbon(alpha=0.2, linetype = "blank") +
  #scale_fill_brewer(palette = "Paired")
p5
ggsave("plots/test20/SBWKilled.png", plot = p5, dpi = 300)



#### By bioclimatic domain ####
SBWKill = data.frame(run=NA, year=NA, phase=NA, cell.id=NA, spp=NA, ny.def=NA, curr.intens.def=NA, bioclim.domain=NA, scn=NA)
#load
for (iscn in c(1,5)){
  SBWKill_scn <- readRDS(paste0("outputs/test7_scn",iscn,"/ap.sbw_results.rds"))$SBWKill
  SBWKill_scn$scn <- paste0("scn_",iscn)
  SBWKill <- rbind(SBWKill, SBWKill_scn)
}
SBWKill <- SBWKill[-1,]

SBWKill$bioclim_merge = ifelse(SBWKill$bioclim.domain<=3, 1, 
                               ifelse(SBWKill$bioclim.domain>=6, 3, 2))
#Taula
taula <- SBWKill %>%
  group_by(scn, run, year, bioclim_merge) %>%
  summarise(km2 = 4*n()) %>%
  group_by(scn, year, bioclim_merge) %>%
  summarise(mean = mean(km2), sd=sd(km2))
taula[is.na(taula)] <- 0

#Plot
#taula = filter(taula, scn=="scn_8")
png(paste0(getwd(), "/plots/test3/SBWKilled_bioclim.png"))
ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("SBWkilled (km2)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank") +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~ bioclim_merge, ncol=3, nrow=1)
dev.off()







### SBWDefol ###############################################################################################
#load
SBWDefol.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA, scn=NA)
for (t in c(20)){
  for (iscn in c(0:9)){
    SBWDefol.sm_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$SBWDefol.sm
    SBWDefol.sm_scn$scn <- paste0("t",t,"_scn",iscn)
    SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)
  }
}
SBWDefol.sm <- SBWDefol.sm[-1,]
SBWDefol.sm = filter(SBWDefol.sm, curr.intens.def!=0)

#Taula
taula <- SBWDefol.sm %>%
  group_by(scn, run, year,curr.intens.def) %>%
  summarise(n_obs = sum(ncell)) %>%
  group_by(scn, year, curr.intens.def) %>%
  summarise(mean = mean(n_obs), sd=sd(n_obs))
taula[is.na(taula)] <- 0
taula$curr.intens.def<-as.factor(taula$curr.intens.def)

#Plot
p6 = ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("SBWDefol (ncell)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")+
  #geom_ribbon(alpha=0.2, linetype = "blank") +
  #scale_fill_brewer(palette = "Paired")+
  facet_wrap(~ curr.intens.def, ncol=3, nrow=1)
p6
ggsave("/plots/test20/SBWDefol_cat.png", plot = p6, dpi = 300)


#### Totes les intensitats juntes ####
taula3 <- SBWDefol.sm %>%
  group_by(scn, run, year) %>%
  summarise(n_obs = sum(ncell)) %>%
  group_by(scn, year) %>%
  summarise(mean = mean(n_obs), sd=sd(n_obs))
taula3[is.na(taula3)] <- 0

#Plot
p7 = ggplot(data=taula3, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn,color=scn)) + 
  xlab("Year") + ylab("SBWDefol (ncell)")+ 
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Paired")#+
  #geom_ribbon(alpha=0.1, linetype = "blank")+
  #scale_fill_brewer(palette = "Paired")
p7
ggsave("/plots/test20/SBWDefol.png", plot = p7, dpi = 300)


#### Continuity with SBWHist ####
load(file=paste0("C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Codi Nuria/Dades SBW/SBW.defol.y.mask.rdata"))
sbw.defol.hist = group_by(sbw.defol.y, year) %>% summarise(km2=sum(area.km))

sbw.defol.proj = filter(taula3, scn=="t20_scn0")
sbw.defol.proj = sbw.defol.proj[,c(2,3)]
sbw.defol.proj$mean = sbw.defol.proj$mean*4
colnames(sbw.defol.proj) = c("year","km2")

sbw.defol.histproj = rbind(sbw.defol.hist,sbw.defol.proj)

#plot
plot(sbw.defol.histproj$year, sbw.defol.histproj$km2, pch=20, col = ifelse(sbw.defol.histproj$year>2020, "orange", "darkred"), ylim=c(0,350000))
lines(sbw.defol.histproj$year, sbw.defol.histproj$km2, col = ifelse(sbw.defol.histproj$year>2020, "orange", "darkred"))


#### By bioclimatic domain ####
SBWDefol = data.frame(run=NA, year=NA, phase=NA, cell.id=NA, spp=NA, curr.intens.def=NA, bioclim.domain=NA, scn=NA)
#load
for (iscn in c(1,5,9)){
  SBWDefol_scn <- readRDS(paste0("outputs/test7_scn",iscn,"/ap.sbw_results.rds"))$SBWDefol
  SBWDefol_scn$scn <- paste0("scn_",iscn)
  SBWDefol <- rbind(SBWDefol, SBWDefol_scn)
}
SBWDefol <- SBWDefol[-1,]

SBWDefol$bioclim_merge = ifelse(SBWDefol$bioclim.domain<=3, 1, 
                               ifelse(SBWDefol$bioclim.domain>=6, 3, 2))
#Taula
taula <- SBWDefol %>%
  group_by(scn, run, year, bioclim_merge) %>%
  summarise(km2 = 4*n()) %>%
  group_by(scn, year, bioclim_merge) %>%
  summarise(mean = mean(km2), sd=sd(km2))
taula[is.na(taula)] <- 0

#Plot
#taula = filter(taula, scn=="scn_8")
png(paste0(getwd(), "/plots/test3/SBWDefoled_bioclim.png"))
ggplot(data=taula, aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  xlab("Year") +  ylab("SBWkilled (km2)") +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank") +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~ bioclim_merge, ncol=3, nrow=1)
dev.off()



### Canvi spp #################################################

SppDist = data.frame(scn=NA, test=NA, year=NA, run=NA, BOJ=NA, EPN=NA, ERS=NA, NonFor=NA, OTH.FEU.N=NA, OTH.FEU.S=NA, OTH.RES.N=NA, OTH.RES.S=NA, PET=NA, SAB=NA)
for (scn in c(0,5,9)){
  for (t in c(20)){
    for (r in c(1,2,3)){
      for (y in seq(2021, 2100, 1)){
        
        land <- readRDS(paste0("outputs/test",t,"_scn",scn,"/landscape_",r,"run_",y,"t.rds"))
        
        SppDist = rbind(SppDist,
                        cbind(data.frame(scn=scn, test=t, year=y, run=r),
                              pivot_wider(as.data.frame.table(table(land$spp)), names_from = Var1, values_from = Freq)))
        
        print( paste0("scn=",scn, "; test=",t, "; run=",r, "; year=",y))
      
      }
    }  
  }
}
SppDist <- SppDist[-1,]



#load
SppDist = data.frame(run=NA, year=NA, BOJ=NA, EPN=NA, ERS=NA, NonFor=NA, OTH.FEU.N=NA, OTH.FEU.S=NA, OTH.RES.N=NA, OTH.RES.S=NA, PET=NA, SAB=NA, scn=NA)
for (iscn in c(0:9,12)){
  SppDist_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SppChange
  SppDist_scn$scn <- paste0("scn_",iscn)
  SppDist <- rbind(SppDist, SppDist_scn)
}
SppDist <- SppDist[-1,]


#Taula
taula <- SppDist %>%
  group_by(scn, year) %>% #group_by(scn, test, year) 
  summarise(BOJ=mean(BOJ), EPN=mean(EPN), ERS=mean(ERS), NonFor=mean(NonFor), OTH.FEU.N=mean(OTH.FEU.N), OTH.FEU.S=mean(OTH.FEU.S),
            OTH.RES.N=mean(OTH.RES.N), OTH.RES.S=mean(OTH.RES.S), PET=mean(PET), SAB=mean(SAB))
taula[is.na(taula)] <- 0

taula2 <- taula %>%
  pivot_longer(cols = c(BOJ, EPN, ERS, NonFor, OTH.FEU.N, OTH.FEU.S, OTH.RES.N, OTH.RES.S, PET, SAB),
               names_to = "spp", values_to = "nobs")
#taula2$scn = as.character(as.numeric(taula2$scn))

#Plot
p8 = ggplot(taula2, aes(x=year, y=nobs, color=scn)) + 
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Paired")+
  facet_wrap(~ spp, ncol=5, nrow=2)
p8
ggsave("/plots/test20/SppChange.png", plot = p8, dpi = 300)




