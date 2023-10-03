##############################################################################################-
#### Script for analyzing ap-sbw model results ###############################################-
#### 

#library
library(tidyr); library(raster); library(sp); library(ggplot2); library(ggpubr); library(dplyr); library(RColorBrewer); library(lme4); library(car); library(ez)


#Dades sbw
dades.intensity <- read.csv("analisis/TempSeries_SBWDefol_t20.csv")
dades.sbw <- read.csv("analisis/TempSeries_SBW_t20.csv")


#Dades spp
SppChange.sm = data.frame(run=NA, year=NA,BOJ=NA,EPN=NA,ERS=NA,NonFor=NA,OTH.FEU.N=NA,OTH.FEU.S=NA,OTH.RES.N=NA,OTH.RES.S=NA,PET=NA,SAB=NA, scn=NA)
for (iscn in c(0,8,9,12)){
  SppChange_scn = readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SppChange
  SppChange_scn$scn = paste0("scn_",iscn)
  SppChange.sm = rbind(SppChange.sm, SppChange_scn)
  
  print(iscn)
}
SppChange.sm <- SppChange.sm[-1,]

#Dades merge
dades.sbw = right_join(dades.sbw, SppChange.sm, by=c("scn","run","year"))
dades.sbw$propSABEPN = dades.sbw$SAB/dades.sbw$EPN
dades.intensity = right_join(dades.intensity, SppChange.sm, by=c("scn","run","year"))
dades.intensity = filter(dades.intensity, curr.intens.def!=0)



#MODEL SBW-INTENSITY
df_sum <- dades.intensity %>%
  group_by(scn, run, year) %>%
  summarise(sbw_int1 = sum(SBWDefol_ncells[curr.intens.def == 1]),
            sbw_int2 = sum(SBWDefol_ncells[curr.intens.def == 2]),
            sbw_int3 = sum(SBWDefol_ncells[curr.intens.def == 3]),
            SAB = max(SAB)) %>%
  mutate(prop_int3 = sbw_int3 / (sbw_int1 + sbw_int2 + sbw_int3))

anova_resultat <- aov(prop_int3 ~ SAB+scn*year + Error(run/year), data = df_sum)
anova_resultat <- aov(prop_int3 ~ scn*year + Error(run/year), data = df_sum)
summary(anova_resultat)

#plot
ggplot(dades.intensity, aes(x = scn, y = SBWDefol_ncells, fill = factor(curr.intens.def))) +
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "Intensity") +
  labs(x = "Scn", y = "Proportion") +
  theme_classic()

#PostHoc Tukey
anova_resultat <- aov(prop_int3 ~ scn*year + (1|run), data = df_sum)
TukeyHSD(anova_resultat, "scn")
#Plot Tukey
tukey_df <- as.data.frame(TukeyHSD(anova_resultat, "scn")$scn)
tukey_df$comparison <- rownames(tukey_df)

ggplot(tukey_df, aes(x = comparison, y = diff)) +
  geom_bar(stat = "identity", aes(fill = abs(`p adj`) < 0.05)) +
  scale_fill_manual(values = c("grey70", "red")) +
  coord_flip() +
  labs(x = "Comparació de nivells d'scn", y = "Diferència mitjana d'area") +
  theme_bw()



#MODEL SBW-SPP
## SBWDefol
summary(aov(SBWDefol_ncells ~ propSABEPN+scn*year + Error(run/year), data = dades.sbw))
summary(aov(SBWDefol_ncells ~ SAB+scn*year + Error(run/year), data = dades.sbw))
summary(aov(SBWDefol_ncells ~ EPN+scn*year + Error(run/year), data = dades.sbw))

anova_resultat <- aov(SBWDefol_ncells ~ SAB+scn*year + Error(run/year), data = filter(dades.sbw, scn%in%c("scn_8","scn_9")))
summary(anova_resultat)
#Plot
ggplot(dades.sbw, aes(x = year, y = SAB)) + 
  geom_point() +   facet_wrap(~scn) +   labs(x = "SAB", y = "SBW")
ggplot(dades.sbw, aes(x = year, y = SBWDefol_ncells)) + 
  geom_point() +   facet_wrap(~scn) +   labs(x = "SAB", y = "SBW")


#Totes les espècies
mod0 = aov(SBWDefol_ncells ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_0"))
summary(mod0)
mod8 <- aov(SBWDefol_ncells ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_8"))
summary(mod8)
mod9 <- aov(SBWDefol_ncells ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_9"))
summary(mod9)
mod12 <- aov(SBWDefol_ncells ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_12"))
summary(mod12)

mod <- aov(SBWDefol_ncells ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+scn*year + (1|run), data = dades.sbw)
summary(mod)
TukeyHSD(mod, "scn")

#plot
ggplot(dades.sbw, aes(x = SAB, y = SBWDefol_ncells)) + 
  geom_point() + 
  facet_wrap(~scn, nrow = 2) + 
  labs(x = "Àrea de l'espècie 1", y = "SBWDefol_ncells") +
  theme_bw()



## SBWKilled
summary(aov(SBWKill_area ~ propSABEPN+scn*year + Error(run/year), data = dades.sbw))
summary(aov(SBWKill_area ~ SAB+scn*year + Error(run/year), data = dades.sbw))
summary(aov(SBWKill_area ~ EPN+scn*year + Error(run/year), data = dades.sbw))

anova_resultat <- aov(SBWKill_area ~ SAB+scn*year + Error(run/year), data = filter(dades.sbw, scn%in%c("scn_8","scn_9")))
summary(anova_resultat)

#Plot
ggplot(dades.sbw, aes(x = year, y = SBWKill_area)) + 
  geom_point() +   facet_wrap(~scn) +   labs(x = "SAB", y = "SBW")

#Totes les espècies
mod0 = aov(SBWKill_area ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_0"))
summary(mod0)
mod8 <- aov(SBWKill_area ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_8"))
summary(mod8)
mod9 <- aov(SBWKill_area ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_9"))
summary(mod9)
mod12 <- aov(SBWKill_area ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+year + Error(run/year), data = filter(dades.sbw, scn=="scn_12"))
summary(mod12)

mod <- aov(SBWKill_area ~ EPN+ERS+OTH.FEU.N+OTH.FEU.S+OTH.RES.N+PET+SAB+scn*year + (1|run), data = dades.sbw)
summary(mod)
TukeyHSD(mod, "scn")

#plot
ggplot(dades.sbw, aes(x = SAB, y = SBWKill_area)) + 
  geom_point() + 
  facet_wrap(~scn, nrow = 2) + 
  labs(x = "Àrea de l'espècie 1", y = "SBWKill_area") +
  theme_bw()






### DADES - AGE
#dades
dades.sbw <- read.csv("analisis/TempSeries_SBW_t20.csv")
#Dades age
AgeSpp.sm = data.frame(run=NA, year=NA,mgmt.unit=NA,spp=NA,age.class=NA,area=NA, scn=NA)
for (iscn in c(0,8,9,12)){
  AgeSpp_scn = readRDS(paste0("outputs/test20_scn",iscn,"/ap.sbw_results.rds"))$SppByAgeClass
  AgeSpp_scn$scn = paste0("scn_",iscn)
  AgeSpp.sm = rbind(AgeSpp.sm, AgeSpp_scn)
  
  print(iscn)
}
AgeSpp.sm <- AgeSpp.sm[-1,]

#Dades merge
dades.age = right_join(dades.sbw, AgeSpp.sm, by=c("scn","run","year"))

#Arrenge1
dades.age$agespp = paste(dades.age$spp,"_",dades.age$age.class)
dades.age2 <- dades.age %>% spread(key = agespp, value = area, sep = "_")


#Arrange
dades.age$age = ifelse(dades.age$age.class=="C10",10,
                       ifelse(dades.age$age.class=="C30",30,
                              ifelse(dades.age$age.class=="C50",50,
                                     ifelse(dades.age$age.class=="C70",70,
                                            ifelse(dades.age$age.class=="C90",90,
                                                   ifelse(dades.age$age.class=="OLD",110,0))))))


#Model tot junt
##Defoliation
mod <- aov(SBWDefol_ncells ~ age.class*spp + scn + year + (1|run), data = dades.age)
summary(mod)

#SABEPN
mod <- aov(SBWDefol_ncells ~ age*spp + scn + year + (1|run), data = filter(dades.age, spp %in% c("SAB","EPN")))
summary(mod)
plot(age)

ggplot(data = filter(dades.age,spp %in% c("SAB","EPN")), aes(x = age, y = SBWDefol_ncells, color = spp)) +
  geom_point() +
  labs(x = "Age", y = "SBWDefol") +
  scale_color_manual(values = c("red", "blue")) 


