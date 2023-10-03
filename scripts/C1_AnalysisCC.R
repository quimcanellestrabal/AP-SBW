##############################################################################################-
#### Script for analyzing ap-sbw model results ###############################################-
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(ggpubr); library(dplyr); library(RColorBrewer); library(lme4); library(car); library(ez)


##################################### -
#### A1: ANALYSIS CLIMATE CHANGE ####
##### Crear taula ####
berta <- data.frame(cell.id=NA, x=NA, y=NA, year=NA, scn=NA, run=NA)

for (r in c(1:10)){
  for (t in 2021:2100){
    land1 = readRDS(paste0("outputs/test20_scn0_noCC/landscape_",r,"run_",t,"t.rds")) %>% filter(curr.intens.def>0) %>% dplyr::select(c("cell.id","x","y")) #referencia a SBWDef
    #land1 = readRDS(paste0("outputs/test20_scn0_noCC/landscape_",r,"run_",t,"t.rds")) %>% filter(tssbw==0) %>% dplyr::select(c("cell.id","x","y")) #referencia a SBWKill
    if(nrow(land1)>0){land1$year = t}
    if(nrow(land1)>0){land1$scn = "NoCC"}
    if(nrow(land1)>0){land1$run = r}
    land2 = readRDS(paste0("outputs/test20_scn0/landscape_",r,"run_",t,"t.rds")) %>% filter(curr.intens.def>0) %>% dplyr::select(c("cell.id","x","y")) #referencia a SBWDef
    #land2 = readRDS(paste0("outputs/test20_scn0/landscape_",r,"run_",t,"t.rds")) %>% filter(tssbw==0) %>% dplyr::select(c("cell.id","x","y")) #referencia a SBWkill
    if(nrow(land2)>0){land2$year = t}
    if(nrow(land2)>0){land2$scn = "scn0_CC"}
    if(nrow(land2)>0){land2$run = r}
    
    if(nrow(land1)>0){berta = rbind(berta,land1)}
    if(nrow(land2)>0){berta = rbind(berta,land2)}
    
    print(paste(t,r))
  }
}

berta <- berta[-1,]
berta$period = ifelse(berta$year<=2040, "peak20s",
                      ifelse(berta$year>2040 & berta$year<=2070, "peak50s",
                             ifelse(berta$year>2070, "peak80s","nopeak")))

write.csv(berta, "analisis/ClimateChange.csv",row.names = FALSE)
berta <- read.csv("analisis/ClimateChange.csv")
##### Linear Model #####
mod1 <- lm(berta$y ~ berta$year * berta$scn)
summary(mod1)

mod2 <- lm(y ~ year, data = subset(berta, scn == "NoCC"))
summary(mod2)

mod3 <- lm(y ~ year, data = subset(berta, scn == "scn0_CC"))
summary(mod3)

#Comparació de les pendents de les regressions
coef2 <- coef(mod2)[2]
coef3 <- coef(mod3)[2]
diff_coef <- coef2 - coef3
#Calcular l'interval de confiança per a la diferència de coeficients
confint(mod2, "year") - confint(mod3, "year")


##### Model mixt #####
mod_mixt <- lmer(y ~ scn*year + (1|run), data=berta)
summary(mod_mixt)
anova(mod_mixt)
car::Anova(mod_mixt, type="III") #el mateix però amb la llibreria car

mod_mixt2 <- lmer(y ~ year + (1|run), data=subset(berta, scn == "NoCC"))
summary(mod_mixt2)

mod_mixt3 <- lmer(y ~ year + (1|run), data=subset(berta, scn == "scn0_CC"))
summary(mod_mixt3)



##### Plots ####
marti = sample_n(berta, 6000)
plot(marti$year, marti$y, pch = 20, col=as.factor(marti$scn))
abline(mod3, col="red")
abline(mod2, col="black")

boxplot(marti$y ~ marti$period*marti$scn)



#################################-
#### A2: CC SERIES TEMPORALS ####

##### Preparar base de dades #####
#SBWDefol
SBWDefol.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA, scn=NA)

  SBWDefol.sm_scn <- readRDS(paste0("outputs/test20_scn0/ap.sbw_results.rds"))$SBWDefol.sm
  SBWDefol.sm_scn$scn <- paste0("scn_0")
  SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)
  
  SBWDefol.sm_scn <- readRDS(paste0("outputs/test20_scn0_noCC/ap.sbw_results.rds"))$SBWDefol.sm
  SBWDefol.sm_scn$scn <- paste0("scn_0_noCC")
  SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)

SBWDefol.sm <- SBWDefol.sm[-1,]

#Taula ini
taula_ini = data.frame(scn = rep(c("scn_0","scn_0_noCC"), each = 3200),
                       run = rep(seq(1,10), each = 320, times = 2),
                       year = rep(seq(2021,2100),each=4, times=20),
                       curr.intens.def = rep(seq(0,3),each=1,times=1600))
#taula_ini = filter(taula_ini, run%in%c(1,2,3))

#Taula SBWDefol_categories
taula_SBWDefol <- SBWDefol.sm %>%
  group_by(scn, run, year,curr.intens.def) %>%
  summarise(SBWDefol_ncells = sum(ncell))

#Taula_SBWDefol final
taula_SBWDefol = left_join(taula_ini, taula_SBWDefol, by=c("scn","run","year","curr.intens.def"))
taula_SBWDefol$scn = as.factor(as.character(taula_SBWDefol$scn))
taula_SBWDefol$run = as.factor(taula_SBWDefol$run)
taula_SBWDefol$curr.intens.def <- as.factor(taula_SBWDefol$curr.intens.def)
taula_SBWDefol[is.na(taula_SBWDefol)] <- 0

#Taula_Defol sense categories
SBWDefol.sm.2 = filter(SBWDefol.sm, curr.intens.def!=0)
taula_SBWDefol2 <- SBWDefol.sm.2 %>%
  group_by(scn, run, year) %>%
  summarise(SBWDefol_ncells = sum(ncell))

#SBWKill
SBWKill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA, scn=NA)

  SBWKill.sm_scn <- readRDS(paste0("outputs/test20_scn0/ap.sbw_results.rds"))$SBWKill.sm
  SBWKill.sm_scn$scn <- paste0("scn_0")
  SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
  
  SBWKill.sm_scn <- readRDS(paste0("outputs/test20_scn0_noCC/ap.sbw_results.rds"))$SBWKill.sm
  SBWKill.sm_scn$scn <- paste0("scn_0_noCC")
  SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
  
SBWKill.sm <- SBWKill.sm[-1,]

#Taula ini
taula_ini = data.frame(scn = rep(c("scn_0","scn_0_noCC"), each = 800),
                       run = rep(seq(1,10), each=80, times = 2),
                       year = rep(seq(2021,2100),each=1, times=20))
#taula_ini = filter(taula_ini, run%in%c(1,2,3))

#Taula SBWKill
taula_SBWKill <- SBWKill.sm %>%
  group_by(scn, run, year) %>%
  summarise(SBWKill_area = sum(area))

#Taula final
dades = left_join(taula_ini, taula_SBWKill, by=c("scn","run","year"))
dades = left_join(dades, taula_SBWDefol2, by=c("scn","run","year"))
dades$scn <- as.factor(as.character(dades$scn))
dades$run <- as.factor(dades$run)
dades[is.na(dades)] <- 0

write.csv(dades, "analisis/CC_TempSeries_SBW_t7.csv",row.names = FALSE)




##### Analysis Anova repeated measures ####
# Carreguem les dades
dades <- read.csv("analisis/CC_TempSeries_SBW_t7.csv")

dades2 <- dades %>%
  group_by(scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

###### SBWKill ####
#Normalitat de les dades
ggqqplot(dades, "SBWKill_area", facet.by = "year")

# Realitzem l'ANOVA per a dades repetides amb tres factors ***no sé que puta merda he de fer!
anova_resultat <- aov(SBWKill_area ~ scn*year + (1|run), data = dades)
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = dades)
summary(anova_resultat)

#Plot 2 a 2
ggplot(data=filter(dades2), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")



