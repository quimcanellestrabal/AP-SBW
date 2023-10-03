### PLOTS FOR MANUSCRIPT ap-SBW ###

library(tidyr); library(raster); library(sp); library(ggplot2); library(dplyr); library(RColorBrewer); library(lme4); library(car); library(ez)

#setwd("C:/Users/Quim/Documents/Quebec/Codi Nuria/AP-SBW2/")

#Color
COL=c(brewer.pal(9, "YlOrRd"))
plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1))
rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)

#Colors ppt
COL=c(rgb(91/255,155/255,213/255),rgb(237/255,125/255,49/255),rgb(255/255,192/255,0/255),rgb(112/255,173/255,71/255),rgb(165/255,165/255,165/255))
plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1))
rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)



##########################-
### STUDY AREA ###########
load("data/mask.rda")
load("data/landscape.rda")
landscape$sppnum = as.numeric(landscape$spp)
table(landscape$spp, landscape$sppnum)

kk=rasterFromXYZ(landscape[,c(2,3,22)])
projection(kk) <- '+proj=lcc +lat_1=46 +lat_2=60 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

plot(kk)

writeRaster(kk, "C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Manuscripts/Figures/land2020_spp2.tif", format = "GTiff", overwrite=T)


######################################-
### HIST-FUTURE CONTINUITY ###########

#Hist
load(file=paste0("C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Codi Nuria/Dades SBW/SBW.defol.y.mask.rdata"))
sbw.defol.hist = group_by(sbw.defol.y, year) %>% summarise(km2=sum(area.km))
sbw.defol.hist$min=sbw.defol.hist$km2
sbw.defol.hist$max=sbw.defol.hist$km2
sbw.defol.hist$time = "hist"

sbw.defol.hist <- sbw.defol.hist[sbw.defol.hist$year != 2021, ]
sbw.defol.hist <- filter(sbw.defol.hist, year>=2005)
sbw.defol.hist$year = as.numeric(as.character(sbw.defol.hist$year))

#Future
SBWDefol <- readRDS(paste0("outputs/test20_scn0/ap.sbw_results.rds"))$SBWDefol.sm
SBWDefol = filter(SBWDefol, curr.intens.def!=0)
sbw.defol.proj <- SBWDefol %>%
  group_by(run, year) %>%
  summarise(n_obs = sum(ncell)) %>%
  group_by(year) %>%
  summarise(km2 = mean(n_obs)*4, min=min(n_obs)*4, max=max(n_obs)*4,)
sbw.defol.proj[is.na(sbw.defol.proj)] <- 0
sbw.defol.proj$time = "proj"

#Tot junt
sbw.defol.histproj = rbind(sbw.defol.hist,sbw.defol.proj)

#Plot
ggplot(data = sbw.defol.histproj, aes(x = year, y = km2)) +
  xlab("Year") + ylab("SBWDefol (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #geom_point(shape=21, color="black", fill="grey80")+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(aes(ymin = min, ymax = max), alpha=0.4, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  geom_line(aes(colour = time),linewidth = 1.2)

#Projectant un run concret
qq <- SBWDefol %>%
  group_by(run, year) %>%
  summarise(km2 = sum(ncell)*4)

ggplot(data = sbw.defol.histproj, aes(x = year, y=km2/10^3)) +
  xlab("Year") + ylab(expression("Area defoliated"~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #geom_point(shape=21, color="black", fill="grey80")+
  geom_ribbon(aes(ymin = min, ymax = max), alpha=0.4, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  geom_line(data=filter(qq, run==6), aes(x = year, y = km2),linewidth = 1.2, color="#D95F02")+ #indicar aquí el run a plotejar
  scale_color_brewer(palette = "Dark2")+
  geom_line(data=sbw.defol.hist, aes(x = year, y = km2),linewidth = 1.2, color="#1B9E77")


ggplot(data = sbw.defol.histproj, aes(x = year, y = (km2/10^3))) +
  xlab("Year") + ylab(expression("Area defoliated"~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_ribbon(aes(ymin = min, ymax = max), alpha=0.4, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  geom_line(data=filter(qq, run==6), aes(x = year, y = km2),linewidth = 1.2, color="#D95F02")+ 
  scale_color_brewer(palette = "Dark2")+
  geom_line(data=sbw.defol.hist, aes(x = year, y = km2),linewidth = 1.2, color="#1B9E77") +
  scale_y_continuous(labels = function(x) x/1000) +  # Dividir la variable y entre 1000
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))  # Ajustar la mida de les etiquetes i de les labels




########################-
### MAPS SBW ###########
load("data/mask.rda")
load("data/landscape.rda")
landscape <- readRDS(paste0("outputs/test20_scn0/landscape_6run_2051t.rds"))
landscape <- readRDS(paste0("outputs/test20_scn0/landscape_6run_2084t.rds"))

kk=rasterFromXYZ(landscape[,c(2,3,19)]) #21 si és per l'any 2020
projection(kk) <- '+proj=lcc +lat_1=46 +lat_2=60 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

plot(kk)

writeRaster(kk, "C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Manuscripts/Figures/SBW_scn0_r6_2020y.tif", format = "GTiff")
writeRaster(kk, "C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Manuscripts/Figures/SBW_scn0_r6_2051y.tif", format = "GTiff", overwrite=T)
writeRaster(kk, "C:/Users/Quim/OneDrive - Universitat de Girona/Quebec/Manuscripts/Figures/SBW_scn0_r6_2084y.tif", format = "GTiff", overwrite=T)


######################################-
### SBW AND CLIMATE CHANGE ###########

COL = brewer.pal(n = 8, name = 'Dark2')

#Plot with CC
land20_sbw <- filter(readRDS("outputs/test20_scn0/landscape_2run_2021t.rds"), curr.intens.def>0)
land50_sbw <- filter(readRDS("outputs/test20_scn0/landscape_2run_2054t.rds"), curr.intens.def>0)
land90_sbw <- filter(readRDS("outputs/test20_scn0/landscape_2run_2089t.rds"), curr.intens.def>0)

plot(density(land90_sbw$y), main="Scn CC", xlab="latitude", xlim=c(180000,1000000),ylim=c(0,0.000006), lwd=2.0)
polygon(density(land90_sbw$y), col=adjustcolor(COL[1], alpha.f=0.7))
lines(density(land50_sbw$y),lwd=2.0)
polygon(density(land50_sbw$y), col=adjustcolor(COL[2], alpha.f=0.7))
lines(density(land20_sbw$y),lwd=2.0)
polygon(density(land20_sbw$y), col=adjustcolor(COL[3], alpha.f=0.7))
legend("topright", legend=c("2020","2050","2080"), fill=adjustcolor(COL[c(3,2,1)],alpha.f = 0.7), cex=1)
abline(v=mean(land90_sbw$y), col=adjustcolor(COL[1], alpha.f=1), lwd=2, lty=2)
abline(v=mean(land50_sbw$y), col=adjustcolor(COL[2], alpha.f=1), lwd=2, lty=2)
abline(v=mean(land20_sbw$y), col=adjustcolor(COL[3], alpha.f=1), lwd=2, lty=2)

#Plot no CC
land20_sbw <- filter(readRDS("outputs/test20_scn0_noCC/landscape_3run_2021t.rds"),curr.intens.def>0)
land50_sbw <- filter(readRDS("outputs/test20_scn0_noCC/landscape_3run_2053t.rds"),curr.intens.def>0)
land90_sbw <- filter(readRDS("outputs/test20_scn0_noCC/landscape_3run_2084t.rds"),curr.intens.def>0)

plot(density(land90_sbw$y), main="Scn BAU", xlab="latitude", xlim=c(180000,1000000),ylim=c(0,0.000006), lwd=2.0)
polygon(density(land90_sbw$y), col=adjustcolor(COL[1], alpha.f=0.7))
lines(density(land50_sbw$y), lwd=2.0)
polygon(density(land50_sbw$y), col=adjustcolor(COL[2], alpha.f=0.7))
lines(density(land20_sbw$y),lwd=2.0)
polygon(density(land20_sbw$y), col=adjustcolor(COL[3], alpha.f=0.7))
legend("topright", legend=c("2020","2050","2080"), fill=adjustcolor(COL[c(3,2,1)],alpha.f = 0.7),cex=1)
abline(v=mean(land90_sbw$y), col=adjustcolor(COL[1], alpha.f=1), lwd=2, lty=2)
abline(v=mean(land50_sbw$y), col=adjustcolor(COL[2], alpha.f=1), lwd=2, lty=2)
abline(v=mean(land20_sbw$y), col=adjustcolor(COL[3], alpha.f=1), lwd=2, lty=2)




################################-
### SERIES TEMPORALS ###########

#### SBWDefol ###################################################################################
#load
SBWDefol.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA, scn=NA)

for (iscn in c(0:9,12)){
  SBWDefol.sm_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SBWDefol.sm
  SBWDefol.sm_scn$scn <- paste0("S",iscn)
  SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)
  print(iscn)
}
SBWDefol.sm_scn <- readRDS(paste0("outputs/test20_scn0_noCC/ap.sbw_results.rds"))$SBWDefol.sm
SBWDefol.sm_scn$scn <- paste0("BAU")
SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)

SBWDefol.sm <- SBWDefol.sm[-1,]
SBWDefol.sm = filter(SBWDefol.sm, curr.intens.def!=0)

#Taula
taula <- SBWDefol.sm %>%
  group_by(scn, run, year) %>%
  summarise(n_obs = sum(ncell)) %>%
  group_by(scn, year) %>%
  summarise(mean = mean(n_obs*4), sd=sd(n_obs*4), 
            confint=qt(0.975, 10-1) * (sd(n_obs)/sqrt(10))) #canvia 10 per el número de rèpliques
taula[is.na(taula)] <- 0
#taula$scn = ifelse(taula$scn=="S12","S1",taula$scn)

#png("C:/Users/Quim/Documents/Quebec/Manuscripts/Figures/SBWDefol.png"))
ggplot(data=filter(taula, scn%in%c("scn_0","scn_6","scn_7","scn_8","scn_9","scn_10")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("SBWDefol (ncell)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  #theme(plot.background = element_rect(fill = "grey90"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.1, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")
#dev.off()


#Figura 5
p1 = ggplot(data=filter(taula, scn%in%c("S0","BAU")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c(brewer.pal(4,"Dark2")[4], brewer.pal(3,"Dark2")[1]))+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  #scale_fill_brewer(values = c(brewer.pal(4,"Dark2")[4], brewer.pal(3,"Dark2")[1]))+
  theme(legend.position = "none")+
  ggtitle("BAU vs CC") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))

p2 = ggplot(data=filter(taula, scn%in%c("S0","S2","S3")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.25%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 12))

p3 = ggplot(data=filter(taula, scn%in%c("S0","S4","S5")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.5%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 12))

p4 = ggplot(data=filter(taula, scn%in%c("S0","S6","S7")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.75%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 12))

p5 = ggplot(data=filter(taula, scn%in%c("S0","S8","S9")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 1%") + theme(plot.title = element_text(size = 15))+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title.x = element_blank())

p100 = ggplot(data=filter(taula, scn%in%c("S0","S1","S9")), aes(x=year, y=mean/10^3, ymin=(mean-confint)/10^3, ymax=(mean+confint)/10^3, fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area defoliated (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c(palette("Dark2")[1], palette("Dark2")[6], palette("Dark2")[3]))+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("No harvesting vs Harvesting vs Harvesting & Artificial planting (BS)") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 12))



#### SBWKilled ###############################################################################################
SBWKill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA, scn=NA)
#load
for (iscn in c(0:9,12)){
    SBWKill.sm_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SBWKill.sm
    SBWKill.sm_scn$scn <- paste0("S",iscn)
    SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
    print(iscn)
}
SBWKill.sm_scn <- readRDS(paste0("outputs/test20_scn0_noCC/ap.sbw_results.rds"))$SBWKill.sm
SBWKill.sm_scn$scn <- paste0("BAU")
SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)

SBWKill.sm <- SBWKill.sm[-1,]

#Taula
taula <- SBWKill.sm %>%
  group_by(scn, run, year) %>%
  summarise(km2 = sum(area)) %>%
  group_by(scn, year) %>%
  summarise(mean = mean(km2), sd=sd(km2), 
            confint=qt(0.975, 10-1) * (sd(km2)/sqrt(10))) #canvia 10 per el número de rèpliques
taula[is.na(taula)] <- 0
#taula$scn = ifelse(taula$scn=="S12","S1",taula2$scn)

#png("C:/Users/Quim/Documents/Quebec/Manuscripts/Figures/SBWKill.png"))
ggplot(data=filter(taula, scn%in%c("scn0","scn1","scn2","scn3","scn4","scn5")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("SBWKill (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  #theme(plot.background = element_rect(fill = "grey90"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.1, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")
#dev.off()

#Figura 5
p6 = ggplot(data=filter(taula, scn%in%c("S0","BAU")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c(brewer.pal(4,"Dark2")[4], brewer.pal(3,"Dark2")[1]))+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  #scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("BAU vs CC") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))

p7 = ggplot(data=filter(taula, scn%in%c("S0","S2","S3")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.25%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))

p8 = ggplot(data=filter(taula, scn%in%c("S0","S4","S5")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.5%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))

p9 = ggplot(data=filter(taula, scn%in%c("S0","S6","S7")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 0.75%") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))

p10 = ggplot(data=filter(taula, scn%in%c("S0","S8","S9")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("Harvesting rate = 1%") + theme(plot.title = element_text(size = 15))+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title.x = element_blank())

p101 = ggplot(data=filter(taula, scn%in%c("S0","S1","S9")), aes(x=year, y=mean, ymin=(mean-confint), ymax=(mean+confint), fill=scn,color=scn)) + 
  xlab("Year") + ylab("Area killed (km2)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_line(linewidth = 1.2)+
  scale_color_manual(values = c(palette("Dark2")[1], palette("Dark2")[6], palette("Dark2")[3]))+
  geom_ribbon(alpha=0.2, linetype = "blank")+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.position = "none")+
  ggtitle("No harvesting vs Harvesting vs Harvesting & Artificial planting (BS)") + theme(plot.title = element_text(size = 15))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_text(size = 12))




##Figura5
library(gridExtra); library(gridtext)
grid.arrange(p1, p100, p2, p3, p4, p5, ncol = 1)
grid.arrange(p6, p101, p7, p8, p9, p10, ncol = 1)





#### LandChange ###############################################################################################
LandChange <- data.frame(run=NA, year=NA, trans=NA, n_obs=NA, scn=NA)

#load
for (iscn in c(0,8,9)){
  LandChange_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$LandChange.sm
  LandChange_scn$scn <- paste0("scn_",iscn)
  LandChange <- rbind(LandChange, LandChange_scn)
}
LandChange <- LandChange[-1,]

#Taula
taula <- LandChange %>%
  group_by(scn, year, trans) %>%
  summarise(mean = mean(n_obs), sd=sd(n_obs),
            confint=qt(0.975, 10-1) * (sd(n_obs)/sqrt(10))) #canvia 10 per el número de rèpliques)

#Plot
#taula = filter(taula, scn=="scn_8")
#png("C:/Users/Quim/Documents/Quebec/Manuscripts/Figures/Transition.png"))
ggplot(data=filter(taula, trans=="Harv_AP"|trans=="Harv_ForTrans"|trans=="NatTrans"), 
       aes(x=year, y=mean, ymin=(mean-sd), ymax=(mean+sd), fill=scn, color=scn)) + 
  facet_wrap(~ trans, ncol=3, nrow=1)+
  xlab("Year") +  ylab("Transition (ncell)") +
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =element_rect(fill="grey20"))+
  theme(strip.text = element_text(colour = 'white'))+
  geom_line(linewidth = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  geom_ribbon(alpha=0.2, linetype="blank") +  
  scale_fill_brewer(palette = "Dark2")
#dev.off()




#############################-
### BOXPLOT - MANAGEMENT ####

#load
dades <- read.csv("analisis/TempSeries_SBW_t20.csv")

dades$scn <- gsub("scn_", "S", dades$scn)

qq = dades %>% group_by(scn,run) %>% summarize(SBWKill=sum(SBWKill_area), SBWDefol=sum(SBWDefol_ncells*4))
kk = data.frame(scn=c("S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S12"),
                artpl=c("C","C","A","B","A","B","A","B","A","B","C"),
                AP=c("A0","A0","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS","A0"),
                harv=c(0,0.25,0.25,0.5,0.5,0.25, 0.75,0.75,1,1,1),
                harv2=c("0","0.25","0.25","0.5","0.5","0.25", "0.75","0.75","1","1","1")) #trampeta amb s3 i s5
dades = merge(qq,kk, by="scn")

dades = filter(dades, scn!="S1")

#Plot
col=brewer.pal(n = 8, name = "Dark2")

p1 = ggplot(dades, aes(x=harv2, y=(SBWDefol*4/10^6), fill=AP)) + 
  geom_boxplot(position = position_dodge(preserve = "single")) + 
  xlab("Harvesting rate (%)") + ylab(expression("Area"~(10^6~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("A0" = col[1], "A50_BS" = col[3], "A50_TA" = col[2]))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  #geom_hline(yintercept = 5.1, linetype = "dashed", color = "darkred", size = 1.5)+
  guides(fill = "none")

p2 = ggplot(dades, aes(x=harv2, y=(SBWKill/10^3), fill=AP)) + 
  geom_boxplot(position = position_dodge(preserve = "single"))+
  xlab("Harvesting rate (%)") + ylab(expression("Area"~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("A0" = col[1], "A50_BS" = col[3], "A50_TA" = col[2]))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  #geom_hline(yintercept = 55, linetype = "dashed", color = "darkred", size = 1.5)+
  guides(fill = "none")



#Proportion
#load
dada <- read.csv("analisis/TempSeries_SBW_t20.csv")
land1 <- read.csv("proves/sbw_potential.csv")
land1$scn <- paste0("scn_", as.character(land1$scn))
land1 = filter(land1, scn!="scn_1")
land1$scn = ifelse(land1$scn=="scn_12", "scn_1",land1$scn)

d10=merge(dada, land1, by=c("scn","run","year"))

d10$SBWKill_area = d10$SBWKill_area/4
d10$propDefolTOTMat = 100*d10$SBWDefol_ncells/d10$p1
d10$propDefolTOT = 100*d10$SBWDefol_ncells/d10$p4
d10$propKillTOTMat = 100*d10$SBWKill_area/d10$p1
d10$propKillTOT = 100*d10$SBWKill_area/d10$p4

#Acumulat
qq = d10 %>% group_by(scn,run) %>% 
  summarize(propDefolTOTMat=mean(propDefolTOTMat), propDefolTOT=mean(propDefolTOT),
            propKillTOTMat=mean(propKillTOTMat), propKillTOT=mean(propKillTOT))
kk = data.frame(scn=c("scn_0","scn_1","scn_2","scn_3","scn_4","scn_5","scn_6","scn_7","scn_8","scn_9"),
                artpl=c("C","C","A","B","A","B","A","B","A","B"),
                AP=c("A0","A0","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS"),
                harv=c(0,1,0.25,0.5,0.5,0.25, 0.75,0.75,1,1),
                harv2=c("0","1","0.25","0.5","0.5","0.25", "0.75","0.75","1","1")) #trampeta amb s3 i s5
d4 = merge(qq,kk, by="scn")
d4 = filter(d4, propDefolTOTMat<400)

#Plot
p3 = ggplot(d4, aes(x=harv2, y=(propDefolTOT), fill=AP)) + 
  geom_boxplot(position = position_dodge(preserve = "single"))+
  xlab("Harvesting rate (%)") + ylab("Proportion (%)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("A0" = col[1], "A50_BS" = col[3], "A50_TA" = col[2]))+
  #geom_hline(yintercept = 5.9, linetype = "dashed", color = "darkred", size = 1.5)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  guides(fill = "none")

p4 = ggplot(d4, aes(x=harv2, y=(propKillTOT), fill=AP)) + 
  geom_boxplot(position = position_dodge(preserve = "single"))+
  xlab("Harvesting rate (%)") + ylab("Proportion (%)")+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("A0" = col[1], "A50_BS" = col[3], "A50_TA" = col[2]))+
  #geom_hline(yintercept = 0.24, linetype = "dashed", color = "darkred", size = 1.5)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  guides(fill = "none")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)




###################-
### CHANGE SPP ####
#load
SppDist = data.frame(run=NA, year=NA, BOJ=NA, EPN=NA, ERS=NA, NonFor=NA, OTH.FEU.N=NA, OTH.FEU.S=NA, OTH.RES.N=NA, OTH.RES.S=NA, PET=NA, SAB=NA, scn=NA)
for (iscn in c(0:9,12)){
  SppDist_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SppChange
  SppDist_scn$scn <- paste0("scn_",iscn)
  SppDist <- rbind(SppDist, SppDist_scn)
  print(iscn)
}
SppDist_scn <- readRDS(paste0("outputs/test20_scn0_noCC/ap.sbw_results.rds"))$SppChange
SppDist_scn$scn <- paste0("scn_0_noCC")
SppDist <- rbind(SppDist, SppDist_scn)

SppDist <- SppDist[-1,]


#Taula
taula <- SppDist %>%
  group_by(scn, year) %>% #group_by(scn, test, year) 
  summarise('Yellow birch'=mean(BOJ), 'Black spruce'=mean(EPN), 'Sugar maple'=mean(ERS), NonFor=mean(NonFor), OTH.FEU.N=mean(OTH.FEU.N), OTH.FEU.S=mean(OTH.FEU.S),
            OTH.RES.N=mean(OTH.RES.N), OTH.RES.S=mean(OTH.RES.S), 'Trembling aspen'=mean(PET), 'Balsam fir'=mean(SAB))
taula[is.na(taula)] <- 0

taula5 <- taula %>%
  mutate('Other boreal' = OTH.RES.N + OTH.FEU.N, 'Other temperate' = OTH.RES.S + OTH.FEU.S) %>% 
  select(!c(NonFor,OTH.RES.N,OTH.FEU.N,OTH.RES.S,OTH.FEU.S))
  

taula2 <- taula5 %>%
  pivot_longer(cols = c('Yellow birch', 'Black spruce', 'Sugar maple', 'Trembling aspen', 'Balsam fir', 'Other boreal', 'Other temperate'),
               names_to = "spp", values_to = "nobs")
#taula2$scn = as.character(as.numeric(taula2$scn))

#Arrenge taula
taula2$scn <- gsub("scn_", "S", taula2$scn)

taula2 = filter(taula2, scn!="S1")
taula2$scn = ifelse(taula2$scn=="S12","S1",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S0_noCC","SBAU",taula2$scn)

#Arrenge new scenario names
taula2$scn = ifelse(taula2$scn=="S0","CC_H0",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S1","CC_H1_A0",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S2","CC_H0.25_A50_TA",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S3","CC_H0.25_A50_BS",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S4","CC_H0.5_A50_TA",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S5","CC_H0.5_A50_BS",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S6","CC_H0.75_A50_TA",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S7","CC_H0.75_A50_BS",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S8","CC_H1_A50_TA",taula2$scn)
taula2$scn = ifelse(taula2$scn=="S9","CC_H1_A50_BS",taula2$scn)
taula2$scn = ifelse(taula2$scn=="SBAU","BAU_H0",taula2$scn)


#Color
COL=c(brewer.pal(10, "Paired"))
COL2=c("#E7298A",COL[c(1,2,5,6,7,8)],"#1B9E77","#E6AB02",COL[c(9,10)])
plot(NULL, xlim=c(0,length(COL2)), ylim=c(0,1))
rect(0:(length(COL2)-1), 0, 1:length(COL2), 1, col=COL2)

#Plot
p8=ggplot(taula2, aes(x=year, y=nobs*4/10^3, color=scn)) + 
  geom_line(linewidth = 1) +
  xlab("Year") + ylab(expression("Area "~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = COL2) +
  facet_wrap(~ spp, ncol=4, nrow=2, scales = "free_x")+
  theme(strip.text = element_text(size = 12),strip.background = element_rect(fill = "grey65"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  scale_x_continuous(breaks = c(2020, 2100), limits = c(2020, 2100))

p8
ggsave("/plots/test20/SppChange.png", plot = p8, dpi = 300)


#Només escenaris significatius
taula3 = filter(taula2, scn%in%c("CC_H0_A0","CC_H1_A0","CC_H1_A50_TA","CC_H1_A50_BS","BAU_H0_A0"))
COL3=c("#E7298A","#1B9E77","#E6AB02",COL[c(9,10)])

ggplot(taula3, aes(x=year, y=nobs*4/10^3, color=scn)) + 
  geom_line(linewidth = 1) +
  xlab("Year") + ylab(expression("Area "~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = COL3) +
  facet_wrap(~ spp, ncol=4, nrow=2, scales = "free_x")+
  theme(strip.text = element_text(size = 12),strip.background = element_rect(fill = "grey65"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  scale_x_continuous(breaks = c(2020, 2100), limits = c(2020, 2100))

ggplot(filter(taula3, spp%in%c("EPN","PET","SAB")), aes(x=year, y=nobs*4/10^3, color=scn)) + 
  geom_line(linewidth = 1) +
  xlab("Year") + ylab(expression("Area "~(10^3~km^2)))+ 
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(palette("Dark2")[1], palette("Dark2")[6],palette("Dark2")[2], palette("Dark2")[3]))+
  facet_wrap(~ spp, ncol=3, nrow=1, scales = "free_x")+
  theme(strip.text = element_text(size = 12),strip.background = element_rect(fill = "grey65"))+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  scale_x_continuous(breaks = c(2020, 2100), limits = c(2020, 2100))



###########################-
### GRAPHICAL ABSTRACT ####

#load
dades <- read.csv("analisis/TempSeries_SBW_t20.csv")

dades$scn <- gsub("scn_", "S", dades$scn)

qq = dades %>% group_by(scn,run) %>% summarize(SBWKill=mean(SBWKill_area), SBWDefol=mean(SBWDefol_ncells*4))
kk = data.frame(scn=c("S0","S1","S2","S3","S4","S5","S6","S7","S8","S9","S12"),
                artpl=c("C","C","A","B","A","B","A","B","A","B","C"),
                AP=c("A0","A0","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS","A50_TA","A50_BS","A0"),
                harv=c(0,0.25,0.25,0.5,0.5,0.25, 0.75,0.75,1,1,1),
                harv2=c("0","0.25","0.25","0.5","0.5","0.25", "0.75","0.75","1","1","1")) #trampeta amb s3 i s5
dades = merge(qq,kk, by="scn")

dades = filter(dades, scn%in%c("S0","S12","S8","S9"))

dades$scn2 = ifelse(dades$scn=="S0","No management",
                   ifelse(dades$scn=="S12","Harvesting",
                          ifelse(dades$scn=="S8","Artificial Planting Trembling aspen",
                                 ifelse(dades$scn=="S9","Artificial Planting Black spruce",dades$scn))))


#Colors ppt
COL=c(rgb(91/255,155/255,213/255),rgb(237/255,125/255,49/255),rgb(255/255,192/255,0/255),rgb(112/255,173/255,71/255),rgb(165/255,165/255,165/255))
COL <- brewer.pal(4, "Dark2")[c(8, 2, 1, 3)]

plot(NULL, xlim=c(0,length(COL)), ylim=c(0,1))
rect(0:(length(COL)-1), 0, 1:length(COL), 1, col=COL)

#Plot
ggplot(dades, aes(x=scn, y=(SBWDefol*4/10^3), fill=scn)) + 
  geom_boxplot()+
  xlab("") + ylab(expression("Area"~(10^3~km^2)))+ 
  theme(axis.text.x = element_blank(),  axis.ticks.x = element_blank())+
  theme(panel.background=element_rect(fill = "white", colour = "grey10", linewidth=1.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = COL)+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  geom_hline(yintercept = 61, linetype = "dashed", color = "darkred", size = 1.5)+
  theme(legend.position = "bottom", legend.direction = "horizontal")

