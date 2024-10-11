#############################################
### Preparing gene expression data from PSYCHENCDOE BrainSpan Dev
### based on Yash Patel, modified by Zhijie Liao

#############################################
#############################################
#############################################
library(tidyverse)
setwd("~/Desktop/CNV_brain")

#rm(list=ls())
# Display the current working directory
#workingDir = "/Users/yash/Downloads/";
#setwd(workingDir);

library(WGCNA);
options(stringsAsFactors = FALSE);

# data downloaded from  http://development.psychencode.org/#
Data = read_table("data/bulkRNA/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt");
#Data[1:3,1:3]
Data<-as.data.frame(Data)
rownames(Data) <- Data$Geneid
Data <- Data[,-1]

# Take a quick look at what is in the data set:
dim(Data);#60155 genes, 607 sample
names(Data);

# get just the expr and transpose 
datExpr0 = as.data.frame(t(Data[,]));

# Loading clinical trait data
traitData = read.csv("data/bulkRNA/mRNA-seq_Sample_metadata.csv");

dim(traitData)
head(traitData)

mydat <- data.frame(matrix(NA, ncol= 1, nrow= 607))
names(mydat) <- "fullcode"
mydat$fullcode <- rownames(datExpr0)
mydat$Braincode <- gsub("\\..*","",mydat$fullcode)
merdat <- merge(mydat, traitData, by = "Braincode", all = T)

# remove columns that hold information we do not need.
allTraits = merdat
dim(allTraits)
names(allTraits)

datExpr <- datExpr0
datTraits <- allTraits

#############################################
#############################################
#############################################

library(WGCNA)
options(stringsAsFactors = FALSE);

#lnames = load(file = "raw_psychENCODE_data.RData");
#lnames
############################################

# Only get data for samples in neoctx, and 11 ROI fully. 
datTraits$fullcode
repo <- datTraits
repo$roi <- gsub(".*\\.","",datTraits$fullcode)
table(repo$roi)
head(repo)
neoctx <- c("DFC", "VFC", "MFC", "OFC", "M1C", "IPC", "S1C", "ITC", "A1C", "STC", "V1C")
repo.fil <- repo[which(repo$roi %in% neoctx),]
table(repo.fil$Window, repo.fil$roi) #
Nsamp <- table(repo.fil$Window, repo.fil$roi); rownames(Nsamp) <- paste("W", rownames(Nsamp), sep = "")
Nsamp; 
repo.fil <- repo[which(repo$roi %in% neoctx & repo$Window < 5 & repo$Window > 1),];table(repo.fil$Window, repo.fil$roi)

tokeep <- repo.fil$fullcode # keep neoctx and Window < 5 [between the range of scRNAseq dat]

orig.Expr <- datExpr; orig.trait <- datTraits # keep the orig in env 
datExpr <- datExpr[(which(rownames(datExpr) %in% tokeep)),] # filtered out the un-needed data
datTraits <- datTraits[(which(datTraits$fullcode %in% tokeep)),]
dim(datExpr)

# remove non-expressed genes in datExpr 
head(names(datExpr))
sum(rownames(datExpr) == datTraits$fullcode)

meanexpr <- colMeans(datExpr)
length(which(meanexpr > 0.5))
ind.to.keep <- which(meanexpr > 0.5)
datExpr.fil <- datExpr[, ind.to.keep] # filtered by keeps that have mean exp > 0.5 RPKM 
datExpr.fil[1:4,1:5]; dim(datExpr.fil)

## have a check on nonexpressed genes
nonexp <- which(meanexpr <= 0.5)
nonexp.mat<-datExpr[, nonexp]
nonexp.mat[1:10,1:4]

table(nonexp.mat<0.5)

geneID<-data.frame(Geneid= colnames(nonexp.mat))%>%separate(col = Geneid, into = c("ENSG", "gene"), remove = F)
head(geneID)
filter(geneID, gene=="GJA8")

# transform log2 expr data
for ( i in 1:ncol(datExpr.fil)) {
  datExpr.fil[,i]   <- log2(datExpr.fil[,i]+1)
}
datExpr.fil[1:4,1:4]
### remove covariate effects using empiralbayesLM function 
names(datTraits) #pH not correctd for since there are some NAs 
datTraits$ROI <- gsub("*.*\\.", "" , datTraits$fullcode) # need to keep ROI 
datTraits$SubjID <- gsub("\\....", "" , datTraits$fullcode) # need to keep ROI 
head(datTraits)
table(is.na(datTraits$pH))
#save(datTraits, datExpr.fil, file = "data/bulkRNA/filtered_nonEXP_prenatal14Donor_April20_2023.RData")
###
## COVARIATE adjustment 
# Bayes method 
names(datTraits)
str(datTraits)
covar.toremove.noday <- datTraits[, c("Sex", "Hemisphere", "RIN","Ethnicity", "Sequencing.Site")] # no age/days regressed out 

covar.tokeep <- as.data.frame(datTraits[,c("ROI")])
rownames(covar.tokeep) <- rownames(datTraits) <- datTraits$fullcode

# run bayeslm function 

datExpr.adj <- empiricalBayesLM(data = datExpr.fil, removedCovariates = covar.toremove.noday, retainedCovariates = covar.tokeep, robustPriors = F, verbose = 2)
adjExp.noday <- as.data.frame(datExpr.adj$adjustedData)

# gene ids from raw data is extracted for changing the names in the next step 
ensembl <- gsub("*\\|.*", "" , (names(datExpr.fil)))
symb <- gsub(".*\\|", "" , (names(datExpr.fil)))
tormuniq <- which(duplicated(symb)) # unique gene symbols 

datExpr.fil <- datExpr.fil[,-tormuniq]
names(datExpr.fil) <- gsub(".*\\|", "" , (names(datExpr.fil)))

adjExp.noday <- adjExp.noday[,-tormuniq]
names(adjExp.noday) <- gsub(".*\\|", "" , (names(adjExp.noday)))

rownames(adjExp.noday) == datTraits$fullcode # double check 
datExpr.fil$ROI <- datTraits$ROI
adjExp.noday$ROI <- datTraits$ROI

# getting the median expression per ROI 
median.expr <- datExpr.fil[,] %>% group_by(ROI) %>% summarise(across(everything(), median))
median.expr.adj.noday <- adjExp.noday[,] %>% group_by(ROI) %>% summarise(across(everything(), median))

mean.expr <- datExpr.fil[,] %>% group_by(ROI) %>% summarise(across(everything(), mean))
mean.expr.adj.noday <- adjExp.noday[,] %>% group_by(ROI) %>% summarise(across(everything(), mean))

mean(datExpr.fil$HOPX[which(datExpr.fil$ROI == "DFC")]) == mean.expr$HOPX[mean.expr$ROI == "DFC"]
median(datExpr.fil$HOPX[which(datExpr.fil$ROI == "DFC")]) == median.expr$HOPX[median.expr$ROI == "DFC"]



