#' ---
#' title: "OLIG2 Analysis"
#' author: "Dr. Stefano Berto"
#' output:
#'  html_document:
#'    theme: united
#'    highlight: tango
#' ---

#' We used HTseq-count to quantify the RNA-seq libraries against GRCh37.87 (protein coding only).
#' We computed normalised exon-level counts (CPM), removing the low abundance reads using a condition cutoff. 
#' Only autosomal genes are retained.
#' In addition to the biological and technical covariates we calculated Surrogate Variables (SVs) based on SVA. 
#' We then used limma to detect schizophrenia differential gene expression.

suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DMwR))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(d3heatmap))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(GGally))
source("Utility_Functions.R")

#' load the data
load("INPUT_DATA/INPUT_OLIG_EXON_SZproj.RData")

# Suppress warnings 
options(warn = -1)
# Create output directory. 	
suffix <- "OLIG2"
output_folder <- paste("OUTPUT_", suffix, "/", sep = "")
dir.create(output_folder,recursive=TRUE,showWarnings = FALSE)

#' Filter the data
location=read.table("INPUT_DATA/Homo_sapiens.GRCh37.87_Protein_Coding_GeneBody.bed")
filtBed=location[grep("chrY|chrX",location$V1),]
data <- log2CPM_OLIG[!(rownames(log2CPM_OLIG) %in% unique(filtBed$V4)),]
class <- unique(TRAITS_OLIG[,1])
filter=apply(data, 1, function(x) (all(x[grep(paste(class[1]),names(x))] > 0) | all(x[grep(paste(class[2]),names(x))] > 0)))
dat <- data[filter,]
pd = TRAITS_OLIG

#' Mutate covariate into factorial
pd$Diagnosis = as.factor(pd$Diagnosis)
pd$Sex = as.factor(pd$Sex)
pd$BrainBank = as.factor(pd$BrainBank)
pd$Hemis = as.factor(pd$Hemis)
pd$Age = as.factor(pd$Age)
pd$Pmi = as.factor(pd$Pmi)

#' Quantile normalization
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
write.table(p, file = paste(output_folder,suffix, "_CPM_QuantNorm", ".txt",sep = ""),sep="\t",quote=F)

#' Variance Explained
var <- VarExp(p,pd,5,FALSE)
pdf("OUTPUT_OLIG2/Variance_Explained_OLIG2.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

#' Register cluster
cl <- makeCluster(3)
registerDoParallel(cl)

#' Define the model
#' Age/Pmi/RIN are continuous so model it as a fixed effect
#' Diagnosis/Sex/Hemis/BrainBank are both categorical, so model them as random effects
form <- ~ (1|Diagnosis) + (1|Sex) + (1|BrainBank) + RIN + (1|Hemis) + (1|Age) + (1|Pmi)
varPart <- fitExtractVarPartModel(p, form, pd)
vp <- varPart[order(varPart$Diagnosis, decreasing=TRUE),]

#' Plot Variance
```{r, message = FALSE}
plotPercentBars(vp[1:50,]) + 
theme_classic()
```

#' PCA
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Diagnosis=pd$Diagnosis)
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)

#' Plot PCA 
```{r, message = FALSE}
ggscatter(PCi, x = "PC1", y = "PC2",color = "Diagnosis",palette=c("steelblue","orange"), shape = 21, size = 3)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() +
ylim(-100,100)+
xlim(-100,100)
```

#' Create model matrix
mod <- model.matrix(~Diagnosis+Sex+BrainBank+Hemis+Age+Pmi+RIN, data =pd)
mod0 <- model.matrix(~Sex+BrainBank+Hemis+Age+Pmi+RIN, data = pd)

#' Create SVA object with 100 permuation 
svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100,method="two-step")
pdSv <- cbind(mod,svaobj$sv)

#' Combine the original pheno data with SVs
pdSVs<- cbind(pd,svaobj$sv)

#' Check all covariates
```{r, message = FALSE}
ggpairs(pdSVs, mapping = aes(color = Diagnosis)) + 
theme_bw()
```

#' Clean with the LM regression only Covariate
pd_sva=cbind(pd[c(-1)],svaobj$sv)
pd_sva <- droplevels(pd_sva)
avebeta.lm<-lapply(1:nrow(p), function(x){
  lm(unlist(p[x,])~.,data = pd_sva)
})
residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
pAdj<-residuals+matrix(apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(pAdj)<-rownames(p)
rownames(residuals) <- rownames(p)
write.table(pAdj, file = paste(output_folder,suffix, "_cleanp_lm_sv", ".txt",sep = ""),sep="\t",quote=F)
write.table(residuals, file = paste(output_folder,suffix, "_residuals_lm_sv", ".txt",sep = ""),sep="\t",quote=F)

#' Robust FIT with limma
fitLM = lmFit(p,pdSv,method="robust");
fitEb = eBayes(fitLM);
OLIG2_DGE = topTable(fitEb, coef = "DiagnosisSchizo",number=nrow(p));
write.table(OLIG2_DGE, file = paste(output_folder,suffix, "_EXON_LIMMA_DGE", ".txt",sep = ""),sep="\t",quote=F)

#' DGE table
DT::datatable(OLIG2_DGE, options = list(pageLength = 10))

#' Save into a .RData
OLIG2_EXON_DATA=list(OLIG2_ExpQuant=p,OLIG2_ExpCleaned_lm=pAdj,OLIG2_pData=pdSVs,OLIG2_Dge_All=OLIG2_DGE)
save(OLIG2_EXON_DATA,file="OUTPUT_OLIG2/OLIG2_EXON_DATA.RData")

#' To run the CrossValidation please do: R CMD BATCH --vanilla CrossValidation_OLIG2.R

#' sessionInfo
sessionInfo()

