##############################################################################################-
#### Script for analyzing ap-sbw model results ###############################################-
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(ggpubr); library(dplyr); library(RColorBrewer); library(lme4); library(car); library(ez)


##############################-
#### A5: SERIES TEMPORALS ####

##### Preparar base de dades #####
#SBWDefol
SBWDefol.sm = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA, t_scn=NA)
for (t in c(7,8)){
  for (iscn in c(1:10)){
    SBWDefol.sm_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$SBWDefol.sm
    SBWDefol.sm_scn$t_scn <- paste0("t",t,"_scn",iscn)
    SBWDefol.sm <- rbind(SBWDefol.sm, SBWDefol.sm_scn)
    
    print(paste(t,iscn))
  }
}
SBWDefol.sm <- SBWDefol.sm[-1,]

#Taula ini
taula_ini = data.frame(test = rep(paste0("t", c(4,6)), each=32000),
                       scn = rep(paste0("scn", seq(1, 10)), each = 3200, times=2),
                       run = rep(seq(1,10), each = 320, times = 20),
                       year = rep(seq(2021,2100),each=4, times=200),
                       curr.intens.def = rep(seq(0,3),each=1,times=16000))
taula_ini$t_scn <- paste0(taula_ini$test, "_", taula_ini$scn)
#taula_ini = filter(taula_ini, run%in%c(1,2,3))

#Taula SBWDefol_categories
taula_SBWDefol <- SBWDefol.sm %>%
  group_by(t_scn, run, year,curr.intens.def) %>%
  summarise(SBWDefol_ncells = sum(ncell))

#Taula_SBWDefol final
taula_SBWDefol = left_join(taula_ini, taula_SBWDefol, by=c("t_scn","run","year","curr.intens.def"))
taula_SBWDefol$scn = as.factor(as.character(taula_SBWDefol$scn))
taula_SBWDefol$run = as.factor(taula_SBWDefol$run)
taula_SBWDefol$curr.intens.def <- as.factor(taula_SBWDefol$curr.intens.def)
taula_SBWDefol[is.na(taula_SBWDefol)] <- 0


#Taula_Defol sense categories
SBWDefol.sm.2 = filter(SBWDefol.sm, curr.intens.def!=0)
taula_SBWDefol2 <- SBWDefol.sm.2 %>%
  group_by(t_scn, run, year) %>%
  summarise(SBWDefol_ncells = sum(ncell))



#SBWKill
SBWKill.sm = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA, t_scn=NA)
for (t in c(7,8)){
  for (iscn in c(1:10)){
    SBWKill.sm_scn <- readRDS(paste0("outputs/test",t,"_scn",iscn,"/ap.sbw_results.rds"))$SBWKill.sm
    SBWKill.sm_scn$t_scn <- paste0("t",t,"_scn",iscn)
    SBWKill.sm <- rbind(SBWKill.sm, SBWKill.sm_scn)
    
    print(paste(t,iscn))
  }
}

SBWKill.sm <- SBWKill.sm[-1,]

#Taula ini
taula_ini = data.frame(test = rep(paste0("t", c(7,8)), each=8000),
                       scn = rep(paste0("scn", seq(1, 10)), each = 800, times=2),
                       run = rep(seq(1,10), each = 80, times = 20),
                       year = rep(seq(2021,2100),each=1, times=200))
taula_ini$t_scn <- paste0(taula_ini$test, "_", taula_ini$scn)
#taula_ini = filter(taula_ini, run%in%c(1,2,3))

#Taula SBWKill
taula_SBWKill <- SBWKill.sm %>%
  group_by(t_scn, run, year) %>%
  summarise(SBWKill_area = sum(area))

#Taula final
dades = left_join(taula_ini, taula_SBWKill, by=c("t_scn","run","year"))
dades = left_join(dades, taula_SBWDefol2, by=c("t_scn","run","year"))
dades$scn <- as.factor(as.character(dades$scn))
dades$run <- as.factor(dades$run)
dades[is.na(dades)] <- 0

write.csv(dades, "analisis/PREM_TempSeries_SBW_t7.csv",row.names = FALSE)



##### Analysis Anova repeated measures ####
# Carreguem les dades
dades <- read.csv("analisis/PREM_TempSeries_SBW_t7.csv")
#dades = filter(dades, scn=="scn_1"|scn=="scn_2"|scn=="scn_5"|scn=="scn_6"|scn=="scn_9"|scn=="scn_10"|scn=="scn_12"|scn=="scn_13")


###### SBWKill ####

## SBWKill a Test8
anova_resultat <- aov(SBWKill_area ~ scn*year + (1|run), data = filter(dades, test=="t8"))
summary(anova_resultat)
#Plot
dades2 <- dades %>%
  filter(test=="t8") %>%
  group_by(scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

ggplot(data=filter(dades2), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")

#PostHoc Tukey
TukeyHSD(anova_resultat, "scn")
#Plot Tukey
tukey_df <- as.data.frame(TukeyHSD(anova_resultat, "scn")$scn)
tukey_df$comparison <- rownames(tukey_df)

ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_bar(stat = "identity", aes(fill = abs(`p adj`) < 0.02)) +
  scale_fill_manual(values = c("grey70", "red")) +
  coord_flip() +
  labs(x = "Comparació de nivells d'scn", y = "Diferència mitjana d'area") +
  theme_bw()


## Test7 vs Test8
anova_resultat <- aov(SBWKill_area ~ test*scn*year + (1|run), data = dades)
summary(anova_resultat)

## Comparació 2 a 2
#t7-t8
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn1"|t_scn=="t8_scn1")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn2"|t_scn=="t8_scn2")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn3"|t_scn=="t8_scn3")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn4"|t_scn=="t8_scn4")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn5"|t_scn=="t8_scn5")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn6"|t_scn=="t8_scn6")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn7"|t_scn=="t8_scn7")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn8"|t_scn=="t8_scn8")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn9"|t_scn=="t8_scn9")))
summary(aov(SBWKill_area ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn10"|t_scn=="t8_scn10")))


#plot
dades3 <- dades %>%
  group_by(t_scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(t_scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

ggplot(data=filter(dades3,t_scn=="t7_scn1"|t_scn=="t8_scn1"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn2"|t_scn=="t8_scn2"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn3"|t_scn=="t8_scn3"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn4"|t_scn=="t8_scn4"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn5"|t_scn=="t8_scn5"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn6"|t_scn=="t8_scn6"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn7"|t_scn=="t8_scn7"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn8"|t_scn=="t8_scn8"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn9"|t_scn=="t8_scn9"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")
ggplot(data=filter(dades3,t_scn=="t7_scn10"|t_scn=="t8_scn10"), aes(x=year, y=SBWKill_mean, ymin=(SBWKill_mean-SBWKill_sd), ymax=(SBWKill_mean+SBWKill_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWkilled (km2)")




###### SBWDefol ####

## SBWDefol a Test8
anova_resultat <- aov(SBWDefol_ncells ~ scn*year + (1|run), data = filter(dades, test=="t8"))
summary(anova_resultat)
#Plot
dades2 <- dades %>%
  filter(test=="t8") %>%
  group_by(scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

ggplot(data=filter(dades2), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=scn, color=scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")

#PostHoc Tukey
TukeyHSD(anova_resultat, "scn")
#Plot Tukey
tukey_df <- as.data.frame(TukeyHSD(anova_resultat, "scn")$scn)
tukey_df$comparison <- rownames(tukey_df)

ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_bar(stat = "identity", aes(fill = abs(`p adj`) < 0.02)) +
  scale_fill_manual(values = c("grey70", "red")) +
  coord_flip() +
  labs(x = "Comparació de nivells d'scn", y = "Diferència mitjana d'area") +
  theme_bw()


## Test7 vs Test8
anova_resultat <- aov(SBWDefol_ncells ~ test*scn*year + (1|run), data = dades)
summary(anova_resultat)

## Comparació 2 a 2
#t7-t8
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn1"|t_scn=="t8_scn1")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn2"|t_scn=="t8_scn2")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn3"|t_scn=="t8_scn3")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn4"|t_scn=="t8_scn4")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn5"|t_scn=="t8_scn5")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn6"|t_scn=="t8_scn6")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn7"|t_scn=="t8_scn7")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn8"|t_scn=="t8_scn8")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn9"|t_scn=="t8_scn9")))
summary(aov(SBWDefol_ncells ~ t_scn*year + (1|run), data = filter(dades, t_scn=="t7_scn10"|t_scn=="t8_scn10")))


#plot
dades3 <- dades %>%
  group_by(t_scn, run, year) %>%
  summarise(SBWKill_area = sum(SBWKill_area), SBWDefol_ncells = sum(SBWDefol_ncells)) %>%
  group_by(t_scn, year) %>%
  summarise(SBWKill_mean = mean(SBWKill_area), SBWKill_sd=sd(SBWKill_area), SBWDefol_mean=mean(SBWDefol_ncells), SBWDefol_sd=sd(SBWDefol_ncells))

ggplot(data=filter(dades3,t_scn=="t7_scn1"|t_scn=="t8_scn1"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn2"|t_scn=="t8_scn2"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn3"|t_scn=="t8_scn3"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn4"|t_scn=="t8_scn4"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn5"|t_scn=="t8_scn5"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn6"|t_scn=="t8_scn6"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn7"|t_scn=="t8_scn7"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn8"|t_scn=="t8_scn8"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn9"|t_scn=="t8_scn9"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")
ggplot(data=filter(dades3,t_scn=="t7_scn10"|t_scn=="t8_scn10"), aes(x=year, y=SBWDefol_mean, ymin=(SBWDefol_mean-SBWDefol_sd), ymax=(SBWDefol_mean+SBWDefol_sd), fill=t_scn, color=t_scn)) + 
  geom_line(linewidth = 1.2) + geom_ribbon(alpha=0.1, linetype = "blank") +  xlab("Year") +  ylab("SBWDefoled (ncells)")



