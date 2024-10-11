## run virtual Ontology
## based on Yash Patel, modified by Zhijie Liao
### Sept17-2023, with corrected scores
library(tidyverse)
library(WGCNA)
### CNV effect size profile
setwd("~/Documents/CNV_brain/")
load("output/mega-analysis_sumStats_5scores_Sept17-2023.RData")

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

### median expression data
load("NeuroCHARGE_Virtual_Ontogeny_scripts/Sharing additional gene ontogeny scripts/median_expr_prenatal_11ctxROI_W2-W4.Rdata")
all.genes <- colnames(median.expr.adj.noday[,-1])
#### cell-specific data.
load("NeuroCHARGE_Virtual_Ontogeny_scripts/Sharing additional gene ontogeny scripts/Bhaduri2020_expression.ready.fin.Rdata")
head(exp_lvl5)
exp_lvl5 <- exp_lvl5 %>% dplyr::filter(Lvl4!="Outlier") %>% dplyr::filter(Lvl4!="Red_blood_cells") %>%
  mutate(Lvl4= gsub("_Neuron", "", Lvl4))

### derive all CNVeffect-profile ~ gene-exp correlation
cnv.exp.cor<-function(cnv.df, cnv, exp.df, cell.df){
  dfcor <- data.frame(matrix(NA, nrow = length(all.genes), ncol = 4),stringsAsFactors=F) 
  names(dfcor) <- c("DX", "Gene", "PrenatalR", "PrenatalPval")
  dfcor[,"DX"] <- cnv
  dfcor[,"Gene" ] <- all.genes
  roi.df<-dplyr::select(cnv.df, "regions")# make all ROIs in the EXP and CNV datasets are in order and the same.
  exp.df<-left_join(roi.df, exp.df, by=c("regions"="ROI"))
  print(exp.df$regions==cnv.df$regions)
  cor <- bicorAndPvalue(x = exp.df[-1], y = cnv.df[,cnv], maxPOutliers = 0.1) 
  dfcor[,3] <- as.vector(cor$bicor)
  dfcor[,4] <- as.vector(cor$p)
  hist(cor$bicor)
  x <- inner_join(dfcor, cell.df, by = "Gene")
  print(head(x))
  x
}

##lof_inv
del_lof_inv.exp<-cnv.exp.cor(del.surf.mega, "Sum_loeuf_inv_DEL", median.expr.adj.noday, exp_lvl5)

dup_lof_inv.exp<-cnv.exp.cor(dup.surf.mega, "Sum_loeuf_inv_DUP", median.expr.adj.noday, exp_lvl5)

###### Permutation test to derive p-values, and 95%CI of null hypothesis
get.pval<-function(all.cor, topN, bootN){
  alldat <- all.cor %>% group_by(Lvl4) %>% arrange(desc(specificity)) %>% dplyr::slice(1:topN)#only topN genes
  #genes.out.top.spec <- setdiff(unique(all.cor$Gene), alldat$Gene) ## topN genes
  out <- alldat%>% summarize(mean_cor = mean(PrenatalR), mean_cor_abs = mean(abs(PrenatalR)))# mean of cor for each celltype
  print(out)
  ## calculate permutated cor:
  boot <- list()
  for ( i in 1:bootN){
    boot[[i]] <- all.cor %>% group_by(Lvl4) %>% arrange(desc(specificity)) %>% dplyr::slice(-(1:topN)) %>%##sample from (12215-200) genes
      dplyr::slice_sample(n=topN, replace = F) %>% summarise(mean_boot_cor = mean(PrenatalR)) %>% ungroup()
  }
  boot.merge <- dplyr::bind_rows(boot); 
  boot.ci <- boot.merge %>% group_by(Lvl4) %>% 
    summarize(blci = quantile(mean_boot_cor, 0.025),
              buci = quantile(mean_boot_cor, 0.975),
              mean_under_null=mean(mean_boot_cor))
  pval.boot <- inner_join(out, boot.merge, by = c("Lvl4"))
  pval.boot <- pval.boot %>% group_by(Lvl4) %>% mutate(mean_under_null=mean(mean_boot_cor)) %>% ungroup()
  pval.boot$mean_cor_raw <- pval.boot$mean_cor
  pval.boot$mean_cor <- pval.boot$mean_cor_raw - pval.boot$mean_under_null
  pval.boot$mean_boot_cor <- pval.boot$mean_boot_cor - pval.boot$mean_under_null
  pval <- pval.boot %>% group_by(Lvl4) %>% summarise(
    pval = sum(abs(c(unique(mean_cor), mean_boot_cor)) >= abs(unique(mean_cor)))/(bootN+1)) %>% ungroup()
  boot.out<-inner_join(out, boot.ci, by = c("Lvl4") )%>%
    inner_join(pval, by ="Lvl4")%>%
    mutate(fdr = p.adjust(pval, method = "fdr"))
  print(boot.out)
  alldat<-left_join(alldat, boot.out, by="Lvl4")
  viron.out<-list(alldat, boot.out)
}


## deletion- 1/LOEUF
del_lof_inv.vo<-get.pval(all.cor= del_lof_inv.exp, topN=200, bootN=10000)
del_lof_inv.vo[[2]]
#sensitivity analysis, with top 100, and top 300 cell genes
del_lof_inv.vo100<-get.pval(all.cor= del_lof_inv.exp, topN=100, bootN=10000)
del_lof_inv.vo300<-get.pval(all.cor= del_lof_inv.exp, topN=300, bootN=10000)

##duplication- 1/LOEUF
dup_lof_inv.vo<-get.pval(all.cor= dup_lof_inv.exp, topN=200, bootN=10000)
dup_lof_inv.vo[[2]]
#sensitivity analysis, with top 100, and top 300 cell genes
dup_lof_inv.vo100<-get.pval(all.cor= dup_lof_inv.exp, topN=100, bootN=10000)
dup_lof_inv.vo100[[2]]
dup_lof_inv.vo300<-get.pval(all.cor= dup_lof_inv.exp, topN=300, bootN=10000)

