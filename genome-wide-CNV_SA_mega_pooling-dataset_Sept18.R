## use sum of scores after correcting the "na.rm=T" bug
## Sept 16, 2023
#Zhijie
library(tidyverse)
library(psych)
library(broom)
library(broom.mixed)
library(lme4)
library(lmerTest)
setwd("~/Desktop/CNV_brain")

##### surface area,  
load("data/clean_11roi_surfaceArea_allCohorts.RData") 
rois<-unique(sys.11roi$roi)
##SYS
sys.total.area<-sys.11roi%>%select(sub_id, totalArea)%>%distinct()
sys.surf<-pivot_wider(sys.11roi, id_cols = sub_id, names_from = roi, values_from = SA)%>%
  left_join(sys.total.area, by="sub_id")
head(sys.surf)
### SPS
sps.total.area<-sps.11roi%>%select(sub_id, totalArea)%>%distinct()
sps.surf<-pivot_wider(sps.11roi, id_cols = sub_id, names_from = roi, values_from = SA)%>%
  left_join(sps.total.area, by="sub_id")
head(sps.surf)
### IMAGEN
imgen.total.area<-imagen.11roi%>%select(sub_id, totalArea)%>%distinct()
imagen.surf<-pivot_wider(imagen.11roi, id_cols = sub_id, names_from = roi, values_from = SA)%>%
  left_join(imgen.total.area, by="sub_id")
head(imagen.surf)
##UKBB
head(ukbb.11roi)
ukbb.total.area<-ukbb.11roi%>%select(sub_id, totalArea)%>%distinct()
ukbb.surf<-pivot_wider(ukbb.11roi, id_cols = sub_id, names_from = roi, values_from = SA)%>%
  left_join(ukbb.total.area, by="sub_id")
head(ukbb.surf)

## combine
surf<-rbind(mutate(sys.surf, cohort="SYS"),
            mutate(sps.surf, cohort="SPS"),
            mutate(imagen.surf, cohort="imagen"),
            mutate(ukbb.surf, cohort="UKBB")  )
head(surf)
##
load("data/FINAL_INDIVIDUAL_SCORES_nCNV_sumLOEUF_15Sept2023_Zhijie.RData")
load("data/imagen/clean_demographic_IMAGEN.RData")##IMAGEN demograph
add.sub0<-function(df){  ##make sub_id (xxx) into same length as "sub-00000xxx" for IMAGEN
  df$sub_id<-str_pad(df$sub_id,12,side = "left", pad = "0")
  df$sub_id<-paste0("sub-",df$sub_id)
  df
}
##SYS +SPS
sys.cnv<-FINAL_INDIVIDUAL_SCORES.3%>%filter(Cohort=='SYS')%>%
  mutate(sub_id=paste(tolower(FID), IID, sep = "_"))%>%# derive id same as brain data 
  select(sub_id, FID, age_months, sex, PC1:PC10,  count_nCNVs_DEL, count_nCNVs_DUP,
         nGenes_DEL,  nGenes_DUP, Sum_oe_lof_upper_DEL,Sum_oe_lof_upper_DUP, 
         Sum_loeuf_inv_DEL, Sum_loeuf_inv_DUP,Sum_oe_lof_DEL, Sum_oe_lof_DUP)%>%
  filter(age_months<400)%>%mutate(site="sys")
head(sys.cnv)
sps.cnv<-FINAL_INDIVIDUAL_SCORES.3%>%filter(Cohort=='SYS')%>%
  mutate(sub_id=paste(tolower(FID), IID, sep = "_"))%>%# derive id same as brain data 
  select(sub_id, FID, age_months, sex, PC1:PC10,  count_nCNVs_DEL, count_nCNVs_DUP,
         nGenes_DEL,  nGenes_DUP, Sum_oe_lof_upper_DEL,Sum_oe_lof_upper_DUP, 
         Sum_loeuf_inv_DEL, Sum_loeuf_inv_DUP,Sum_oe_lof_DEL, Sum_oe_lof_DUP)%>%
  filter(age_months>400)%>%mutate(site="sps")

###IMAGEN
imagen.cnv<-FINAL_INDIVIDUAL_SCORES.3%>%filter(Cohort=='Imagen')%>%# row 1026 id: 62936394-2-660Wq
  separate(col="individual", into = c("sub_id","teq"))%>%# derive id same as brain data 
  select(sub_id, FID, age_months, sex, PC1:PC10, count_nCNVs_DEL, count_nCNVs_DUP,
         nGenes_DEL,  nGenes_DUP, Sum_oe_lof_upper_DEL,Sum_oe_lof_upper_DUP, 
         Sum_loeuf_inv_DEL, Sum_loeuf_inv_DUP,Sum_oe_lof_DEL, Sum_oe_lof_DUP)%>%add.sub0()%>%
  inner_join(select(demo, sub_id, site), by="sub_id")
head(imagen.cnv)  
##UKBB
ukbb.cnv<-FINAL_INDIVIDUAL_SCORES.3%>%filter(Cohort=='UKBB')%>%
  mutate(sub_id=FID)%>% #FID= ind ID
  select(sub_id, FID, age_months, sex, PC1:PC10, count_nCNVs_DEL, count_nCNVs_DUP,
         nGenes_DEL,  nGenes_DUP, Sum_oe_lof_upper_DEL,Sum_oe_lof_upper_DUP, 
         Sum_loeuf_inv_DEL, Sum_loeuf_inv_DUP,Sum_oe_lof_DEL, Sum_oe_lof_DUP)%>%
  filter(sub_id%in% ukbb.surf$sub_id)%>%mutate(site="ukbb")
head(ukbb.cnv)
colnames(sys.cnv)==colnames(imagen.cnv);colnames(imagen.cnv)==colnames(ukbb.cnv)
cnv<-rbind(sys.cnv,sps.cnv, imagen.cnv, ukbb.cnv)
cnv[is.na(cnv)]<-0
cnv$age_months[cnv$age_months==0]<-NA
head(cnv)
###
surf.cnv<-inner_join(surf, cnv, by="sub_id")
dim(surf.cnv)
colnames(surf.cnv)
carrier_noncarrier<-surf.cnv%>%mutate(del_dup= nGenes_DEL+ nGenes_DUP)
##
surf.cnv%>% filter(nGenes_DEL>0 )%>%group_by(sex)%>%
  summarise(N=n(),  ageMean= mean(age_months/12, na.rm=T),
            ageSD= sd(age_months/12, na.rm=T))

surf.cnv%>% filter(nGenes_DUP>0 )%>%group_by(sex)%>%
  summarise(N=n(),  ageMean= mean(age_months/12, na.rm=T),
            ageSD= sd(age_months/12, na.rm=T))
##plot
## Number of carrier with nGene del/dup
library(ggpubr)
colnames(surf.cnv)
ngene_ind<-surf.cnv%>%select(sub_id, nGenes_DEL, nGenes_DUP)%>%
  pivot_longer(cols = 2:3, names_to = "TYPE", values_to = "Ngene_ind")%>%
  mutate(TYPE=gsub("nGenes_", "", TYPE))%>%filter(Ngene_ind>0)
head(ngene_ind)
loeuf_ind<-surf.cnv%>%select(sub_id, Sum_loeuf_inv_DEL, Sum_loeuf_inv_DUP)%>%
  pivot_longer(cols = 2:3, names_to = "TYPE", values_to = "loeuf")%>%
  mutate(TYPE=gsub("Sum_loeuf_inv_", "", TYPE))%>%filter(loeuf>0)
head(loeuf_ind)
max(loeuf_ind$loeuf)
np1<-ggplot(ngene_ind,aes(x=Ngene_ind, fill=TYPE, color=TYPE, alpha=TYPE))+
  geom_histogram(binwidth = 1, position="dodge")+
  theme_light(base_size = 12)+xlab("Number of genes deleted/duplicated ")+
  ylab("Number of individuals")+
  scale_x_continuous(breaks = seq(0, 60, 5))+
  scale_alpha_manual(values=c(1, 0.5))+
  theme(legend.position = c(0.85, 0.7))+labs(fill=NULL, color=NULL, alpha=NULL) 
np2<-ggplot(loeuf_ind,aes(x=loeuf, fill=TYPE, color=TYPE, alpha=TYPE))+
  geom_histogram(binwidth = 1, position="dodge")+
  theme_light(base_size = 12)+xlab("Sum of 1/LOEUF")+
  ylab("Number of individuals")+
  scale_x_continuous(breaks = seq(0, 80, 5))+
  scale_alpha_manual(values=c(1, 0.5))+
  theme(legend.position = c(0.85, 0.7))+labs(fill=NULL, color=NULL, alpha=NULL) 
tiff("outplot/Number_of_carrier_with_nGene_LOEUF_del-dup.tiff", units = "cm",width = 13, height = 12, res=500)
ggarrange(np2, np1, nrow = 2)
dev.off()


ggplot(surf.cnv%>%filter(nGenes_DEL>0), aes(x=nGenes_DEL))+
  geom_bar()+
  geom_bar(data=surf.cnv%>%filter(nGenes_DUP>0), aes(x=nGenes_DUP), fill="blue3")
ggplot(surf.cnv%>%filter(nGenes_DUP>0), aes(x=nGenes_DUP))+geom_bar()


surf.cnv%>%filter(site=="sps" | site=="ukbb" |site=="sys")%>%
  group_by(site)%>%summarise(N=n(), nF=sum(sex=="F"), nM=sum(sex=="M"), 
                             ageMean= mean(age_months/12, na.rm=T),
                             ageSD= sd(age_months/12, na.rm=T))

surf.cnv%>%filter(site!="sps" & site!="ukbb" &site!="sys")%>%mutate(cohort="IMAGEN")%>%
  group_by(cohort)%>%
  summarise(N=n(), nF=sum(sex=="F"), nM=sum(sex=="M"), 
                             ageMean= mean(age_months/12, na.rm=T),
                             ageSD= sd(age_months/12, na.rm=T))


#### pooling all  #################
## adjust: sex, age, PCs, 
# with family ID as random effect## with family ID as random effect
mix.mod<-function(df, roi.list, cnv.var){## with family ID as random effect
  df<-as.data.frame(df)
  all.stat<-tibble()
  for (cv in cnv.var) {
    sstat<-tibble()
    for (i in roi.list) {
      m<-lmer(scale(df[,i]) ~df[, cv]+sex+site+age_months+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ (1|FID), df)
      s<-tidy(m, effects="fixed")%>%as.data.frame()%>%mutate(regions=i)
      sstat<-rbind(sstat, s)
      print(paste(cv, i))
    }
    sstat<-filter(sstat, term=="df[, cv]")%>%select(-term, -effect)%>%mutate(cnv=cv)
    all.stat<-rbind(all.stat, sstat)
  }
  all.stat<-all.stat%>%mutate(regions=factor(regions, levels = c(roi.order, "totalArea")), N=nrow(df))
  print(head(all.stat))
  all.stat
}
del.var<-c("nGenes_DEL",  "Sum_loeuf_inv_DEL")#"Sum_oe_lof_DEL", 
dup.var<-c("nGenes_DUP",  "Sum_loeuf_inv_DUP")#"Sum_oe_lof_DUP",
roi.order<-c("DFC", "VFC","OFC", "MFC", "ITC","STC", "IPC", "M1C", 'S1C', "A1C",  "V1C")# 

##deletion
del.surf.mix<-mix.mod(surf.cnv, roi.list = c(rois), cnv.var = del.var)%>%mutate(method="mix")%>%#
  mutate(p.level= case_when(p.value<0.05 ~ "p.value<0.05",
                            p.value>=0.05 ~ "non.sig"))
##duplication
dup.surf.mix<-mix.mod(surf.cnv, roi.list = c(rois), cnv.var = dup.var)%>%mutate(method="mix")%>%#
  mutate(p.level= case_when(p.value<0.05 ~ "p.value<0.05",
                            p.value>=0.05 ~ "non.sig"))

