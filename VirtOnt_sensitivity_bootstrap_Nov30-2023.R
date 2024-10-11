
##  run bootstrap for CNV effects ~ celltype gene expression
### Nov 30
### Zhijie

library(tidyverse)
library(WGCNA)
setwd("~/Desktop/CNV_brain")
## phenotype profile
load("output/mega-analysis_sumStats_5scores_Sept17-2023.RData")
head(del.surf.mix); table(del.surf.mix$cnv)
head(dup.surf.mix);
del.surf.mega<-pivot_wider(del.surf.mix, id_cols = "regions",
                           names_from = "cnv", values_from = "estimate")%>%
  dplyr::filter(regions!="totalArea")%>%
  mutate_at(.vars = 2:ncol(.), .funs = abs)##so higher values, bigger effect size
del.surf.mega

dup.surf.mega<-pivot_wider(dup.surf.mix, id_cols = "regions",
                           names_from = "cnv", values_from = "estimate")%>%
  dplyr::filter(regions!="totalArea")%>%
  mutate_at(.vars = 2:ncol(.), .funs = abs)##so higher values, bigger effect size
dup.surf.mega

### clinic CNV effect size profile
load("output/sumamry_stats_clinicCNV_BC+Cardiff+UKBB_sameCtrl_Nov17.RData")

del.cnv<-all.cnv.stat.del%>%separate("cnv", into = c("CNV", "type"), remove = F)%>%
  pivot_wider(id_cols = "regions", names_from = "CNV", values_from = "estimate")%>%
  dplyr::filter(regions!="totalArea")
colnames(del.cnv)<-gsub("1q","q1", colnames(del.cnv)); colnames(del.cnv)<-gsub("16p","p16", colnames(del.cnv))
colnames(del.cnv)<-gsub("22q","q22", colnames(del.cnv))
del.cnv<-del.cnv%>%mutate(q1=-1*q1, q22=-1*q22)##so higher values, bigger negative (positive) effect
del.cnv

## duplication
head(all.cnv.stat.dup)
dup.cnv<-all.cnv.stat.dup%>%separate("cnv", into = c("CNV", "type"), remove = F)%>%
  pivot_wider(id_cols = "regions", names_from = "CNV", values_from = "estimate")%>%
  dplyr::filter(regions!="totalArea")
colnames(dup.cnv)<-gsub("1q","q1", colnames(dup.cnv)); colnames(dup.cnv)<-gsub("16p","p16", colnames(dup.cnv))
colnames(dup.cnv)<-gsub("22q","q22", colnames(dup.cnv))
dup.cnv<-dup.cnv%>%mutate(p16=-1*p16)##so higher values, bigger negative (positive) effect
dup.cnv

##### obtain gene list for each celltype
load("output/lof_inv.vo_effect-1600Cell_gene_cor.RData")
head(lof_inv.vo[[1]])
dim(lof_inv.vo[[1]])
cell.type<-unique(lof_inv.vo[[1]]$Lvl4)

######load expression data and donors infor
load('data/bulkRNA/filtered_nonEXP_adjExp.noday_prenatal14Donor_May10.Rdata')

adjExp.noday[1:3,1:3]
adjExp.noday[1:3, 19700: 19703]

sample.inf<-data.frame(id=rownames(adjExp.noday))%>%separate(col="id", into = c("donor", "ROI"))
head(sample.inf)
donors<-unique(sample.inf$donor)

## exp for bootstrap
adjExp.1600<-adjExp.noday[, c(lof_inv.vo[[1]]$Gene)]%>%as.data.frame()%>%
  rownames_to_column("id")%>%separate(col="id", into = c("donor", "ROI"))
adjExp.1600[1:3,1:3]
dim(adjExp.1600)

get.boot.mean<-function(cnv.df, cnv, cell.types, cellgene.1600.df){
  ##1. resample donor
  boot.donor<-sample(donors, size = length(donors), replace = T)# donors= 14 donors
  boot.donor.df<-tibble() 
  for (bd in boot.donor) {# resample donor
    bd.df<-filter(adjExp.1600,donor==bd) # adjExp.1600: donor_region * 1600 genes
    boot.donor.df<-rbind(boot.donor.df, bd.df)
  }
  median.expr0 <- boot.donor.df[-1] %>% group_by(ROI) %>% summarise(across(everything(), median))%>%ungroup# median of bootstrapped exp
  roi.df<-dplyr::select(median.expr0, "ROI")# make all ROIs in the EXP and CNV datasets are in order and the same.
  cnv.df<-left_join(roi.df, cnv.df, by=c("ROI"="regions"))
  #print(cnv.df$ROI==median.expr$ROI) ## double-check ROI order of exp and phe are the same
  ##for each celltype: resample genes  and then calculate phe ~exp  correlation
  ##2. resample genes for each celltype
  boot.cell.mean<-NULL
  for (coi in cell.type) {
    ## resample gene 
    cell.gene<-filter(cellgene.1600.df, Lvl4==coi) #gene list
    boot.gene<-sample(cell.gene$Gene, size = 200, replace = T)#resample gene for each celltype
    #boot.df<-boot.donor.df[,c("ROI",boot.gene)]
    median.expr<-median.expr0[,c("ROI",boot.gene)]
    #median.expr <- boot.df %>% group_by(ROI) %>% summarise(across(everything(), median))%>%ungroup# median of bootstrapped exp
    ### correlate with phenotype ## checked ROI order of exp and phe are the same in the above
    cor <- bicorAndPvalue(x = median.expr[-1], y = cnv.df[,cnv], maxPOutliers = 0.1) 
    mean.cor<-mean(cor$bicor)
    boot.cell.mean<-c(boot.cell.mean, mean.cor)
  }
  boot.cell.mean.df<-data.frame(cell=cell.types, boot.mean=boot.cell.mean)%>%column_to_rownames("cell")
  boot.cell.mean.df
}

#x<-get.boot.mean(cnv.df=del.surf.mega, cnv= "Sum_loeuf_inv_DEL", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
#get.boot.mean(cnv.df=del.cnv, cnv= "q1", cell.types= cell.type)
#get.boot.mean(cnv.df=dup.cnv, cnv= "q1", cell.types= cell.type)

### get bootstrap distribution with 1000 repeatition for lof_inv
Nboot=1000
#del
Nboot.mean.del<-tibble_row()
for (i in 1:Nboot) {
  set.seed(i*2)
  boot.mean.df<-get.boot.mean(cnv.df=del.surf.mega, cnv= "Sum_loeuf_inv_DEL", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  Nboot.mean.del<-cbind(Nboot.mean.del, boot.mean.df)
  print(paste("boot",i))
}
Nboot.mean.del[1:4,1:4]

#dup
Nboot.mean.dup<-tibble_row()
for (i in 1:Nboot) {
  set.seed(i*2)
  boot.mean.df<-get.boot.mean(cnv.df=dup.surf.mega, cnv= "Sum_loeuf_inv_DUP", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  Nboot.mean.dup<-cbind(Nboot.mean.dup, boot.mean.df)
  print(paste("boot",i))
}
Nboot.mean.dup[1:4,1:4]

## save for later analysis and plot
lof.inv.boot.vo.del<-t(Nboot.mean.del)%>%as.data.frame()
lof.inv.boot.vo.dup<-t(Nboot.mean.dup)%>%as.data.frame()
#save(lof.inv.boot.vo.del, lof.inv.boot.vo.dup, file = "data/VirOnt_DEL_DUP_lof_inv_boot_Nov-2023.RData")

#t(Nboot.mean.del)%>%as.data.frame()%>%describe(., skew = F)%>%as.data.frame()

### get bootstrap distribution with 1000 repeatition for clinic CNV del
Nboot=1000
q1.del.boot.vo<-tibble_row()
p16.del.boot.vo<-tibble_row()
q22.del.boot.vo<-tibble_row()

for (i in 1:Nboot) {
  set.seed(i*2)
  #q1
  boot.mean.df<-get.boot.mean(cnv.df=del.cnv, cnv= "q1", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  q1.del.boot.vo<-cbind(q1.del.boot.vo, boot.mean.df)
  #p16
  boot.mean.df<-get.boot.mean(cnv.df=del.cnv, cnv= "p16", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  p16.del.boot.vo<-cbind(p16.del.boot.vo, boot.mean.df)
  #q22
  boot.mean.df<-get.boot.mean(cnv.df=del.cnv, cnv= "q22", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  q22.del.boot.vo<-cbind(q22.del.boot.vo, boot.mean.df)
  print(paste("boot",i))
}
q1.del.boot.vo<-q1.del.boot.vo%>%t()%>%as.data.frame()
p16.del.boot.vo<-p16.del.boot.vo%>%t()%>%as.data.frame()
q22.del.boot.vo<-q22.del.boot.vo%>%t()%>%as.data.frame()

q1.dup.boot.vo<-tibble_row()
p16.dup.boot.vo<-tibble_row()
q22.dup.boot.vo<-tibble_row()
for (i in 1:Nboot) {
  set.seed(i*2)
  #q1
  boot.mean.df<-get.boot.mean(cnv.df=dup.cnv, cnv= "q1", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  q1.dup.boot.vo<-cbind(q1.dup.boot.vo, boot.mean.df)
  #p16
  boot.mean.df<-get.boot.mean(cnv.df=dup.cnv, cnv= "p16", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  p16.dup.boot.vo<-cbind(p16.dup.boot.vo, boot.mean.df)
  print(paste("boot",i))
  #q22
  boot.mean.df<-get.boot.mean(cnv.df=dup.cnv, cnv= "q22", cell.types= cell.type, cellgene.1600.df=lof_inv.vo[[1]])
  q22.dup.boot.vo<-cbind(q22.dup.boot.vo, boot.mean.df)
  print(paste("boot",i))
}
q1.dup.boot.vo<-q1.dup.boot.vo%>%t()%>%as.data.frame()
p16.dup.boot.vo<-p16.dup.boot.vo%>%t()%>%as.data.frame()
q22.dup.boot.vo<-q22.dup.boot.vo%>%t()%>%as.data.frame()

save(lof.inv.boot.vo.del, lof.inv.boot.vo.dup,
     q1.del.boot.vo, q1.dup.boot.vo,
     p16.del.boot.vo, p16.dup.boot.vo, 
     q22.del.boot.vo, q22.dup.boot.vo,
     file = "data/VirOnt_all_CNV_boot_Nov-2023.RData")

#############3
################## get p value
get.boot.p<-function(x){
  lci<-quantile(x, probs = 0.025)%>%as.numeric()
  hci<-quantile(x, probs = 0.975)%>%as.numeric()
  xm<-mean(x)
  xc<-x-xm
  pval<-(sum(abs(xc) > abs(xm))+1)/(length(x)+1)
  p.df<-data.frame(bootm= xm, lci = lci, hci=hci, pval=pval)
  p.df
}
get.boot.p(lof.inv.boot.vo.del[,"Endothelial"])

#
get.boot.p.cell<-function(boot.df, cnv){
  boot.p<-tibble()
  for (coi in cell.type) {
    boot.p<-rbind(boot.p, get.boot.p(boot.df[, coi]) )
  }
  boot.p<-boot.p%>%mutate(cell=cell.type,
                          pfdr=p.adjust(pval, method = "fdr"),
                          cnv= cnv)
  boot.p
}


#######lof inv
del.lof.boot.p<-get.boot.p.cell(lof.inv.boot.vo.del, "genome-wide CNV Del")
dup.lof.boot.p<-get.boot.p.cell(lof.inv.boot.vo.dup, "genome-wide CNV Dup")
### clinic CNV
q1.del.boot.p<-get.boot.p.cell(q1.del.boot.vo, "1q21.1 Del")
p16.del.boot.p<-get.boot.p.cell(p16.del.boot.vo, "16p11.2 Del")
q22.del.boot.p<-get.boot.p.cell(q22.del.boot.vo, "22q11.2 Del")

q1.dup.boot.p<-get.boot.p.cell(q1.dup.boot.vo, "1q21.1 Dup")
p16.dup.boot.p<-get.boot.p.cell(p16.dup.boot.vo, "16p11.2 Dup")
q22.dup.boot.p<-get.boot.p.cell(q22.dup.boot.vo, "22q11.2 Dup")

all.boot.p<-rbind(del.lof.boot.p, dup.lof.boot.p,
                  q1.del.boot.p, p16.del.boot.p, 
                  q22.del.boot.p, q1.dup.boot.p, 
                  p16.dup.boot.p, q22.dup.boot.p)%>%
  mutate(p.level= case_when(pfdr<0.05 ~ "p<0.05",
                            pfdr>=0.05 ~ "non.sig"))%>%
  mutate(cell= gsub("_Neuron", "", cell))

write.csv(all.boot.p, "output/VirtOnt_bootstrap_result_Nov30.csv", row.names = F)
#######
all.boot.p<-rbind(del.lof.boot.p, dup.lof.boot.p,
                  q1.del.boot.p, p16.del.boot.p, 
                  q22.del.boot.p, q1.dup.boot.p, 
                  p16.dup.boot.p, q22.dup.boot.p)%>%
  mutate(p.level= case_when(pval<0.05 ~ "p<0.05",
                            pval>=0.05 ~ "non.sig"))%>%
  mutate(cell= gsub("_Neuron", "", cell))
####plot
all.boot.plot<-rbind(lof.inv.boot.vo.del%>%mutate(CNV="genome-wide CNV Del"),
                     lof.inv.boot.vo.dup%>%mutate(CNV="genome-wide CNV Dup"),
                     q1.del.boot.vo%>%mutate(CNV= "1q21.1 Del"),
                     p16.del.boot.vo%>%mutate(CNV= "16p11.2 Del"),
                     q22.del.boot.vo%>%mutate(CNV= "22q11.2 Del"),
                     q1.dup.boot.vo%>%mutate(CNV= "1q21.1 Dup"),
                     p16.dup.boot.vo%>%mutate(CNV= "16p11.2 Dup"),
                     q22.dup.boot.vo%>%mutate(CNV= "22q11.2 Dup"))

head(all.boot.plot)
cell.oreder<-c("Radial_Glia", "IPC", "Excitatory", "Inhibitory", "Microglia", "OPC","Endothelial", "Mural")
cnv.oreder<-c("genome-wide CNV Del", "1q21.1 Del", "16p11.2 Del", "22q11.2 Del","genome-wide CNV Dup", "1q21.1 Dup", "16p11.2 Dup", "22q11.2 Dup")
all.boot.plot.long<-all.boot.plot%>%pivot_longer(cols = 1:8, names_to = "cell", values_to = "bicor")%>%
  mutate(cell= gsub("_Neuron", "", cell))%>%
  left_join(select(all.boot.p, cnv, cell, p.level), by=c("CNV"="cnv", "cell"))%>%
  mutate(cell=factor(cell, levels = cell.oreder),
         CNV=factor(CNV, levels = cnv.oreder))
  
head(all.boot.plot.long)

tiff("outplot/VirtOnt_bootstrap_plot_pval_Nov30-2023.tiff", units = "cm",width = 20, height = 11, res=300)
ggplot(all.boot.plot.long, aes(x=CNV, y=bicor, color=cell))+
  geom_hline(yintercept = 0)+
  geom_boxplot(outlier.size = 0)+
  facet_wrap(facets = "cell", nrow = 1)+
  geom_point(aes(x=CNV, y=-0.4, size=p.level), shape= 8, color="black")+
  theme_linedraw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5))+
  scale_size_manual(values =c(NA, 1))+
  xlab("")
dev.off()


