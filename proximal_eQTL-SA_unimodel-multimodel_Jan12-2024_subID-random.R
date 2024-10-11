##prenatal eQTL PGS and surface area with UKBB
## proximal genes & Distal genes
## unimodal vs. multimodal
## model: lmer(SA ~ PGS + region + (1|donor-ID))
### 
### 
library(tidyverse)
library(broom)
library(lme4)
library(lmerTest)
library(broom.mixed)
setwd("~/Desktop/CNV_brain")

##proximal (deleted) gene list
load("data/deleted_proximal-distal-genes_Jan11-2024.RData")
proximal_gene<-individual_gene_proximal_del%>%filter(!is.na(prog_high))%>%
  mutate(eqtl_id=paste("Pt_0.01_",gene_id, sep = ""))

##Surface area ##### surface area,  
load("data/clean_11roi_surfaceArea_allCohorts.RData") 
rois<-unique(sys.11roi$roi)
##UKBB
ukbb.total.area<-ukbb.11roi%>%select(sub_id, totalArea)%>%distinct()
ukbb.surf<-pivot_wider(ukbb.11roi, id_cols = sub_id, names_from = roi, values_from = SA)%>%
  left_join(ukbb.total.area, by="sub_id")
####UKBB, demographic data
load("data/FINAL_INDIVIDUAL_SCORES_nCNV_sumLOEUF_15Sept2023_Zhijie.RData")
ukbb.dem<-FINAL_INDIVIDUAL_SCORES.3%>%filter(Cohort=='UKBB')%>%
  filter(is.na(nGenes_DEL) | nGenes_DEL==0)%>%#exclude CNV carrier with >=1 deleted
  filter(is.na(nGenes_DUP) | nGenes_DUP==0)%>%#exclude CNV carrier with >=1 duplicated
  mutate(sub_id=FID)%>% #FID= ind ID
  select(sub_id, FID, age_months, sex, PC1:PC10)%>%
  filter(sub_id%in% ukbb.surf$sub_id)

#eQTL scores
load("data/Proximal_eQTL_scores_Jan11-2024.RData")
proximal.eqtl.scores<-proximal.eqtl.scores%>%filter(sub_id%in% ukbb.dem$sub_id)
#####
ukbb.p01.pgs<-inner_join(ukbb.dem, proximal.eqtl.scores, by="sub_id")%>%inner_join(ukbb.surf, by="sub_id")
#
all_region<-c("DFC", "VFC", "OFC", "MFC", "ITC", "STC", "IPC","M1C", "S1C", "A1C", "V1C")
colnames(ukbb.p01.pgs)
for (rg in c(all_region)) {## scale the SA for each region first,
  ukbb.p01.pgs[,rg]<-scale(ukbb.p01.pgs[,rg])## so in the mixed effect modal, no scaling for SA & no region as covariate
}
allroi.ukbb.p01.pgs<-ukbb.p01.pgs%>%pivot_longer(cols = c(multim_region, unim_region), names_to = "region", values_to = "allROI")
##
pc.mod.mixed<-function(df, roi.list, cnv.var){
  df<-as.data.frame(df)
  all.stat<-tibble()
  for (cv in cnv.var) {
    sstat<-tibble()
    for (i in roi.list) {
      m<-lmer(df[,i] ~scale(df[, cv])+ scale(age_months)+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ (1|sub_id), df)
      s<-tidy(m, effects="fixed")%>%as.data.frame()%>%mutate(regions=i)
      sstat<-rbind(sstat, s)
      x<-summary(m); #print(paste("DF:",x$df))
    }
    sstat<-filter(sstat, term=="scale(df[, cv])")%>%select(-term, -effect)%>%mutate(cnv=cv)
    all.stat<-rbind(all.stat, sstat)
  }
  all.stat<-all.stat%>%mutate(regions=factor(regions, levels = c("multimod", "unimod", "allROI")), N=nrow(df))
  print(head(all.stat))
  all.stat
}
proximal.pgs.var<-colnames(proximal.eqtl.scores[-1])
#get result
proximal.eqtl.allroi<-pc.mod.mixed(allroi.ukbb.p01.pgs, roi.list = "allROI", cnv.var = proximal.pgs.var)
