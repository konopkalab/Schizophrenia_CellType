#' ---
#' title: "Cell-Type Analysis"
#' author: "Dr. Stefano Berto"
#' output:
#'  html_document:
#'    theme: united
#'    highlight: tango
#'  pdf_document: default
#' ---

#' We used HTseq-count to quantify the RNA-seq libraries against GRCh37.87 (protein coding only).
#' We computed normalised exon-level counts (CPM), removing the low abundance reads using a condition cutoff. 
#' Only autosomal genes are retained.
#' In addition to the biological and technical covariates we calculated Surrogate Variables (SVs) based on SVA. 
#' We then used limma to detect schizophrenia differential gene expression.

suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DMwR))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GGally))
source("Utility_Functions.R")

#' load the data
load("INPUT_DATA/INPUT_NEUN_EXON_SZproj.RData")
load("INPUT_DATA/INPUT_OLIG_EXON_SZproj.RData")

exp <- cbind(log2CPM_NEUN,log2CPM_OLIG)
exp_ctl <- exp[,grep("Control",names(exp))]

# Suppress warnings 
options(warn = -1)
# Create output directory. 	
suffix <- "CELLTYPE"
output_folder <- paste("OUTPUT_", suffix, "/", sep = "")
dir.create(output_folder,recursive=TRUE,showWarnings = FALSE)

#' Filter the data
location=read.table("INPUT_DATA/Homo_sapiens.GRCh37.87_Protein_Coding_GeneBody.bed")
filtBed=location[grep("chrY|chrX",location$V1),]
data <- exp_ctl[!(rownames(exp_ctl) %in% unique(filtBed$V4)),]
TRAITS <- rbind(TRAITS_NEUN,TRAITS_OLIG)
TRAITS_CTL <- TRAITS[TRAITS$Diagnosis == "Control",]
class <- c("NeuN","Olig2")
filter=apply(data, 1, function(x) (all(x[grep(paste(class[1]),names(x))] > 0) | all(x[grep(paste(class[2]),names(x))] > 0)))
dat <- data[filter,]
pd = TRAITS_CTL
pd$Diagnosis <- c(rep("NeuN",27),rep("Olig2",22))
write.table(pd, file = paste(output_folder,suffix, "_pd", ".txt",sep = ""),sep="\t",quote=F)
pd <- read.table("OUTPUT_CELLTYPE/CELLTYPE_pd.txt")

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
write.table(p, file = paste(output_folder,suffix, "_Vst_QuantNorm", ".txt",sep = ""),sep="\t",quote=F)

#' Variance Explained by Covariates
```{r, message = FALSE}
var <- VarExp(p,pd,10,FALSE)
plotVarExp(var,"Variance Explained: NeuN") +
ggsave("OUTPUT_NeuN/NEUN_VarExp.pdf", width=6, height=6)
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
mod <- model.matrix(~Diagnosis+Sex+BrainBank+Hemis+Age+RIN, data =pd)
mod0 <- model.matrix(~Sex+BrainBank+Hemis+Age+RIN, data = pd)

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
CELL_DGE = topTable(fitEb, coef = "DiagnosisOlig2",number=nrow(p));
write.table(CELL_DGE, file = paste(output_folder,suffix, "_LIMMA_DGE", ".txt",sep = ""),sep="\t",quote=F)

#' DGE table
DT::datatable(CELL_DGE, options = list(pageLength = 10))

#' Save into a .RData
CELL_DATA=list(CELL_ExpQuant=p,CELL_ExpCleaned_lm=pAdj,CELL_pData=pdSVs,CELL_Dge_All=CELL_DGE)
save(CELL_DATA,file="OUTPUT_CELLTYPE/CELL_DATA.RData")

#' sessionInfo
sessionInfo()
