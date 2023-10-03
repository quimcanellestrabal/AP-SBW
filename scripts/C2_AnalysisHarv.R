##############################################################################################-
#### Script for analyzing ap-sbw model results ###############################################-
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(ggpubr); library(dplyr); library(RColorBrewer); library(lme4); library(car); library(ez)


##############################-
#### A3: SERIES TEMPORALS ####

##### Preparar base de dades #####
#SBWDefol
SBWDefol.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA, scn=NA)
for (iscn in c(0:9,12)){
  SBWDefol.sm_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SBWDefol.sm
  SBWDefol.sm_scn$scn <- paste0("scn_",iscn)
  SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)
  
  print(iscn)
}
SBWDefol.sm <- SBWDefol.sm[-1,]

#Taula ini
taula_ini = data.frame(scn = rep(paste0("scn_", c(seq(0, 14),100)), each = 3200),
                       run = rep(seq(1,10), each = 320, times = 16),
                       year = rep(seq(2021,2100),each=4, times=160),
                       curr.intens.def = rep(seq(0,3),each=1,times=12800))
taula_ini = filter(taula_ini, run%in%c(1:10), scn%in%SBWDefol.sm$scn)

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

write.csv(taula_SBWDefol, "analisis/TempSeries_SBWDefol_t20.csv",row.names = FALSE)

#Taula_Defol sense categories
SBWDefol.sm.2 = filter(SBWDefol.sm, curr.intens.def!=0)
taula_SBWDefol2 <- SBWDefol.sm.2 %>%
  group_by(scn, run, year) %>%
  summarise(SBWDefol_ncells = sum(ncell))



#SBWKill
SBWKill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA, scn=NA)
for (iscn in c(0:9,12)){
  SBWKill.sm_scn <- readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SBWKill.sm
  SBWKill.sm_scn$scn <- paste0("scn_",iscn)
  SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
  
  print(iscn)
}
SBWKill.sm <- SBWKill.sm[-1,]

#Taula ini
taula_ini = data.frame(scn = rep(paste0("scn_", c(seq(0, 14),100)), each = 800),
                       run = rep(seq(1,10), each=80, times = 16),
                       year = rep(seq(2021,2100),each=1, times=160))
taula_ini = filter(taula_ini, run%in%c(1:10), scn%in%SBWKill.sm$scn)

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

write.csv(dades, "analisis/TempSeries_SBW_t20.csv",row.names = FALSE)






##### Analysis Anova repeated measures ####
# Carreguem les dades
dades <- read.csv("analisis/TempSeries_SBW_t20.csv")
#dades = filter(dades, scn=="scn_1"|scn=="scn_2"|scn=="scn_5"|scn=="scn_6"|scn=="scn_9"|scn=="scn_10"|scn=="scn_12"|scn=="scn_13")

dades2 <- dades %>%
  group_by(scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

###### SBWKill ####
#Normalitat de les dades
dades %>%
  group_by(year) %>%
  shapiro_test(SBWKill_area)
ggqqplot(dades, "SBWKill_area", facet.by = "year")


# Realitzem l'ANOVA per a dades repetides amb tres factors ***no sé que puta merda he de fer!
anova_resultat <- aov(SBWKill_area ~ scn*year + (1|run), data = dades)
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = dades)
summary(anova_resultat)
#Plot
ggplot(data=filter(dades2), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")

#Posthoc 2 a 2
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_0"|scn=="scn_12"))
summary(anova_resultat)
#Plot 2 a 2
ggplot(data=filter(dades2, scn%in%c("scn_0","scn_12")), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")

#PostHoc Tukey
anova_resultat <- aov(SBWKill_area ~ scn*year + (1|run), data = dades)
TukeyHSD(anova_resultat, "scn")
#Plot Tukey
tukey_df <- as.data.frame(TukeyHSD(anova_resultat, "scn")$scn)
tukey_df$comparison <- rownames(tukey_df)
tukey_df2 <- tukey_df[grep("scn_0", tukey_df$comparison), ]

ggplot(tukey_df2, aes(x = comparison, y = diff)) +
  geom_bar(stat = "identity", aes(fill = abs(`p adj`) < 0.05)) +
  scale_fill_manual(values = c("grey70", "red")) +
  coord_flip() +
  labs(x = "Comparació de nivells d'scn", y = "Diferència mitjana d'area") +
  theme_bw()

#plot entre scn diferents
ggplot(data=filter(dades2, scn%in%c("scn_0","scn_9","scn_12")), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")



###### Comparació entre PET-EPN ####
#2-3
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_3"|scn=="scn_2"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_2"|scn=="scn_3"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#4-5
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_5"|scn=="scn_4"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_4"|scn=="scn_5"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#6-7
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_7"|scn=="scn_6"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_6"|scn=="scn_7"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#8-9
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_8"|scn=="scn_9"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_8"|scn=="scn_9"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#9-10
anova_resultat <- aov(SBWKill_area ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_9"|scn=="scn_10"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_9"|scn=="scn_10"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")



###### SBWDefol ####
#Normalitat de les dades
dades %>%
  group_by(year) %>%
  shapiro_test(SBWDefol_ncells)
ggqqplot(dades, "SBWDefol_ncells", facet.by = "year")


# Realitzem l'ANOVA per a dades repetides amb tres factors ***no sé que puta merda he de fer!
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + (1|run), data = dades)
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = dades)
summary(anova_resultat)
#Plot
ggplot(data=filter(dades2), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")

#Posthoc 2 a 2
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_12"|scn=="scn_0"))
summary(anova_resultat)
#Plot 2 a 2
ggplot(data=filter(dades2,scn=="scn_0"|scn=="scn_1"|scn=="scn_2"|scn=="scn_4"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")

#PostHoc Tukey
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + (1|run), data = dades)
TukeyHSD(anova_resultat)
#Plot Tukey
tukey_df <- as.data.frame(TukeyHSD(anova_resultat, "scn")$scn)
tukey_df$comparison <- rownames(tukey_df)
tukey_df2 <- tukey_df[grep("scn_12", tukey_df$comparison), ]

ggplot(tukey_df2, aes(x = comparison, y = diff)) +
  geom_bar(stat = "identity", aes(fill = abs(`p adj`) < 0.05)) +
  scale_fill_manual(values = c("grey70", "red")) +
  coord_flip() +
  labs(x = "Comparació de nivells d'scn", y = "Diferència mitjana d'area") +
  theme_bw()
#Plot entre scn diferents
ggplot(data=filter(dades2,scn=="scn_0"|scn=="scn_12"|scn=="scn_9"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")



###### Comparació entre PET-EPN ####
#2-3
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_3"|scn=="scn_2"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_2"|scn=="scn_3"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")
#4-5
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_5"|scn=="scn_4"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_4"|scn=="scn_5"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")
#6-7
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_7"|scn=="scn_6"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_6"|scn=="scn_7"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")
#8-9
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_9"|scn=="scn_8"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_8"|scn=="scn_9"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")
#9-10
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + Error(run/year), data = filter(dades, scn=="scn_9"|scn=="scn_10"))
summary(anova_resultat)
ggplot(data=filter(dades2,scn=="scn_9"|scn=="scn_10"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWKill_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefol (pixel)")



########################-
#### A4: HARV - SBW ####

qq = dades %>% group_by(scn) %>% summarize(SBWKill=mean(SBWKill_area), SBWDefol=mean(SBWDefol_ncells))
kk = data.frame(scn=c("scn_0","scn_1","scn_2","scn_3","scn_4","scn_5","scn_6","scn_7","scn_8","scn_9"),
                harv=c(0,0.25,0.25,0.25,0.5,0.5, 0.75,0.75,1,1))
kk = merge(qq,kk, by="scn")
kk$x=kk$harv
kk$y=kk$SBWDefol
kk$z=kk$SBWKill

#Model SBWDefol
mod <- glm(kk$SBWDefol~kk$harv, family = poisson(link = "log"))
mod <- lm(SBWDefol~harv + I(harv^2), data = kk)
summary(mod)
#Plot SBWDefol
fit2 <- lm(y~poly(x,2,raw=TRUE), data=kk)
plot(kk$x, kk$y, pch=19, xlab='harvest ratio (%)', ylab='SBWdefol (ncells)')
x_axis <- seq(0, 1, length=16)
lines(x_axis, predict(fit2, data.frame(x=x_axis)), col='red')

#Model SBWKill
mod <- lm(SBWKill~harv + I(harv^2), data = kk)
summary(mod)
#Plot SBWKill
fit2 <- lm(z~poly(x,2,raw=TRUE), data=kk)
plot(kk$x, kk$z, pch=19, xlab='harvest ratio (%)', ylab='SBWKilled (km2)')
x_axis <- seq(0, 1, length=16)
lines(x_axis, predict(fit2, data.frame(x=x_axis)), col='red')


#El mateix ben fet:
kk=dades
kk$harv <- 0
kk$harv[kk$scn == "scn_1" | kk$scn == "scn_2" | kk$scn == "scn_3"] <- 0.25
kk$harv[kk$scn == "scn_4" | kk$scn == "scn_5"] <- 0.5
kk$harv[kk$scn == "scn_6" | kk$scn == "scn_7"] <- 0.75
kk$harv[kk$scn == "scn_8" | kk$scn == "scn_9"] <- 1

mod <- lm(SBWDefol_ncells~poly(harv,2,raw=TRUE), data=kk)
summary(mod)
plot(kk$harv, kk$SBWDefol_ncells, pch=19, xlab='harvest ratio (%)', ylab='SBWdefol (ncells)')
harv <- seq(0, 1, length=16)
lines(harv, predict(mod, data.frame(x=harv)), col='red')

mod <- lm(SBWKill_area~poly(harv,2,raw=TRUE), data=kk)
summary(mod)
plot(kk$harv, kk$SBWKill_area, pch=19, xlab='harvest ratio (%)', ylab='SBWKill (km2)')
harv <- seq(0, 1, length=16)
lines(harv, predict(mod, data.frame(x=harv)), col='red')



########################-
#### A5: CC - noCC ####
# Carreguem les dades
dades_CC <- read.csv("analisis/TempSeries_SBW_t20.csv")
dades_CC$Clim = "CC"
dades_noCC <- read.csv("analisis/TempSeries_SBW_t20_noCC.csv")
dades_noCC$Clim = "noCC"

dades = rbind(dades_CC, dades_noCC)
#dades = filter(dades, scn=="scn_1"|scn=="scn_2"|scn=="scn_5"|scn=="scn_6"|scn=="scn_9"|scn=="scn_10"|scn=="scn_12"|scn=="scn_13")

dades$scn2 = paste0(dades$scn,"_",dades$Clim)

dades2 <- dades %>%
  group_by(scn2, scn, Clim, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))


###### SBWKill ####
# Realitzem l'ANOVA per a dades repetides amb tres factors ***no sé que puta merda he de fer!
anova_resultat <- aov(SBWKill_area ~ scn*year*Clim + Error(run/year), data = dades)
summary(anova_resultat)
#Plot
ggplot(data=filter(dades2), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")

ggplot(data=filter(dades2), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=Clim, color=Clim)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")



###### Comparació entre CC-noCC ####
#0
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_0_CC"|scn2=="scn_0_noCC")))
ggplot(data=filter(dades2,scn2=="scn_0_CC"|scn2=="scn_0_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#1
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_1_CC"|scn2=="scn_1_noCC")))
ggplot(data=filter(dades2,scn2=="scn_1_CC"|scn2=="scn_1_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#2
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_2_CC"|scn2=="scn_2_noCC")))
ggplot(data=filter(dades2,scn2=="scn_2_CC"|scn2=="scn_2_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#3
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_3_CC"|scn2=="scn_3_noCC")))
ggplot(data=filter(dades2,scn2=="scn_3_CC"|scn2=="scn_3_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#4
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_4_CC"|scn2=="scn_4_noCC")))
ggplot(data=filter(dades2,scn2=="scn_4_CC"|scn2=="scn_4_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#5
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_5_CC"|scn2=="scn_5_noCC")))
ggplot(data=filter(dades2,scn2=="scn_5_CC"|scn2=="scn_5_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#6
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_6_CC"|scn2=="scn_6_noCC")))
ggplot(data=filter(dades2,scn2=="scn_6_CC"|scn2=="scn_6_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#7
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_7_CC"|scn2=="scn_7_noCC")))
ggplot(data=filter(dades2,scn2=="scn_7_CC"|scn2=="scn_7_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#8
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_8_CC"|scn2=="scn_8_noCC")))
ggplot(data=filter(dades2,scn2=="scn_8_CC"|scn2=="scn_8_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
#9
summary(aov(SBWKill_area ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_9_CC"|scn2=="scn_9_noCC")))
ggplot(data=filter(dades2,scn2=="scn_9_CC"|scn2=="scn_9_noCC"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")



###### SBWDefol ####
# Realitzem l'ANOVA per a dades repetides amb tres factors ***no sé que puta merda he de fer!
anova_resultat <- aov(SBWDefol_ncells ~ scn*year*Clim + Error(run/year), data = dades)
summary(anova_resultat)
#Plot
ggplot(data=filter(dades2), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")

ggplot(data=filter(dades2), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=Clim, color=Clim)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")


###### Comparació entre CC-noCC ####
#0
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_0_CC"|scn2=="scn_0_noCC")))
ggplot(data=filter(dades2,scn2=="scn_0_CC"|scn2=="scn_0_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#1
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_1_CC"|scn2=="scn_1_noCC")))
ggplot(data=filter(dades2,scn2=="scn_1_CC"|scn2=="scn_1_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#2
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_2_CC"|scn2=="scn_2_noCC")))
ggplot(data=filter(dades2,scn2=="scn_2_CC"|scn2=="scn_2_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#3
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_3_CC"|scn2=="scn_3_noCC")))
ggplot(data=filter(dades2,scn2=="scn_3_CC"|scn2=="scn_3_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#4
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_4_CC"|scn2=="scn_4_noCC")))
ggplot(data=filter(dades2,scn2=="scn_4_CC"|scn2=="scn_4_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#5
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_5_CC"|scn2=="scn_5_noCC")))
ggplot(data=filter(dades2,scn2=="scn_5_CC"|scn2=="scn_5_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#6
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_6_CC"|scn2=="scn_6_noCC")))
ggplot(data=filter(dades2,scn2=="scn_6_CC"|scn2=="scn_6_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#7
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_7_CC"|scn2=="scn_7_noCC")))
ggplot(data=filter(dades2,scn2=="scn_7_CC"|scn2=="scn_7_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#8
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_8_CC"|scn2=="scn_8_noCC")))
ggplot(data=filter(dades2,scn2=="scn_8_CC"|scn2=="scn_8_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")
#9
summary(aov(SBWDefol_ncells ~ scn2*year + Error(run/year), data = filter(dades, scn2=="scn_9_CC"|scn2=="scn_9_noCC")))
ggplot(data=filter(dades2,scn2=="scn_9_CC"|scn2=="scn_9_noCC"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn2, color=scn2)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDef (ncells)")


##### HARV - SBW by PREM, CC, noCC ####

dades_CC <- read.csv("analisis/TempSeries_SBW_t20.csv")
dades_CC$Clim = "CC"
dades_noCC <- read.csv("analisis/TempSeries_SBW_t20_noCC.csv")
dades_noCC$Clim = "noCC"
dades_PREM <- read.csv("analisis/TempSeries_SBW_t20_PREM.csv")
dades_PREM$Clim = "PREM"

dades = rbind(dades_CC, dades_noCC, dades_PREM)

#Taula amb totes
qq_CC = dades_CC %>% group_by(scn,Clim) %>% summarize(SBWKill=mean(SBWKill_area), SBWDefol=mean(SBWDefol_ncells))
qq_noCC = dades_noCC %>% group_by(scn,Clim) %>% summarize(SBWKill=mean(SBWKill_area), SBWDefol=mean(SBWDefol_ncells))
qq_PREM = dades_PREM %>% group_by(scn,Clim) %>% summarize(SBWKill=mean(SBWKill_area), SBWDefol=mean(SBWDefol_ncells))

kk = data.frame(scn=c("scn_0","scn_1","scn_2","scn_3","scn_4","scn_5","scn_6","scn_7","scn_8","scn_9"),
                harv=c(0,0.25,0.25,0.25,0.5,0.5, 0.75,0.75,1,1))

kk_CC = merge(qq_CC,kk, by="scn")
kk_noCC = merge(qq_noCC,kk, by="scn")
kk_PREM = merge(qq_PREM,kk, by="scn")

kk = rbind(kk_CC, kk_noCC, kk_PREM)
kk$x=kk$harv
kk$y=kk$SBWDefol
kk$z=kk$SBWKill

#Model
mod1 <- glm(kk$SBWDefol~kk$harv*kk$Clim)
summary(mod1)

mod2 <- glm(kk$SBWKill~kk$harv*kk$Clim)
summary(mod2)

#Plot
ggplot(kk, aes(x = x, y = z, color = Clim)) +
  geom_point(size=3) +  scale_color_discrete(name = "Clim") +
  geom_smooth(aes(group = Clim, color = Clim), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE) +
  labs(x = "harvest ratio (%)", y = "SBWKill (km2)") +
  theme_classic()



