#' ---
#' title: "BrainSeq Analysis"
#' author: "Dr. Stefano Berto"
#' output:
#'  html_document:
#'    theme: united
#'    highlight: tango
#'  pdf_document: default
#' ---

#' Data was downloaded from http://eqtl.brainseq.org/phase1/.
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
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(SummarizedExperiment))
source("Utility_Functions.R")

# Suppress warnings 
options(warn = -1)
# Create output directory. 	
suffix <- "BrainSeq"
output_folder <- paste("OUTPUT_", suffix, "/", sep = "")
dir.create(output_folder,recursive=TRUE,showWarnings = FALSE)

#load the data
#load("INPUT_DATA/rse_exon_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda")
#exp <- as.data.frame(assays(rse_exon)$counts)
#TRAITS_BrainSeq <- as.data.frame(colData(rse_exon))
#rows <- as.data.frame(rowData(rse_exon))

# Subset based on age > 16 
#TRAITS_BrainSeq <- TRAITS_BrainSeq[TRAITS_BrainSeq$Age > 16,]
#exp <- exp[,colnames(exp) %in% rownames(TRAITS_BrainSeq)]
#TRAITS_BrainSeq <- droplevels(TRAITS_BrainSeq)

# Annotate the data and collapse duplicated genes
#exp$ID <- rows$Geneid
#exp$Symbol <- mapIds(EnsDb.Hsapiens.v86,keys=exp$ID, column="GENENAME",keytype="GENEID",multiVals="first")
#exp$BioType <- mapIds(EnsDb.Hsapiens.v86,keys=exp$ID, column="TXBIOTYPE",keytype="GENEID",multiVals="first")
#exp$Chr <- mapIds(EnsDb.Hsapiens.v86,keys=exp$ID, column="SEQNAME",keytype="GENEID",multiVals="first")
#exp <- exp[exp$BioType == "protein_coding" & !(exp$Chr %in% c("X","Y","MT")),]
#exp$BioType <- NULL
#exp$ID <- NULL
#exp$Chr <- NULL
#exp <- na.omit(exp)
#agg <- aggregate(x = exp[, 1:ncol(exp)-1], by = list(ID = exp$Symbol), FUN = "sum", na.rm = T)
#rownames(agg) <- agg$ID
#agg$ID <- NULL
#exp <- agg
#save(exp,TRAITS_BrainSeq,file="BrainSeq_INPUT.RData")

#' Filter the data
load("INPUT_DATA/BrainSeq_INPUT.RData")
dat <- log2(cpm(exp)+1)
class <- unique(TRAITS_BrainSeq$Dx)
filter=apply(dat, 1, function(x) (all(x[grep(paste(class[1]),names(x))] > 0) | all(x[grep(paste(class[2]),names(x))] > 0)))
dat <- dat[filter,]

#' Fix the demographic
pd <- TRAITS_BrainSeq
pd <- pd[c(8,4,5,6,7,55)]
colnames(pd)[1] <- "Diagnosis"
colnames(pd)[6] <- "Pmi"
pd$Age <- ifelse(pd$Age <= 40,"1", ifelse (pd$Age > 40 & pd$Age <= 60 ,"2", ifelse(pd$Age > 60,"3","3")))
pd$Pmi <- ifelse(pd$Pmi <= 12,"1", ifelse (pd$Pmi > 12 & pd$Pmi <= 24 ,"2", ifelse(pd$Pmi > 24,"3","3")))

#' Mutate covariate into factorial
pd$Diagnosis = as.factor(pd$Diagnosis)
pd$Sex = as.factor(pd$Sex)
pd$Age = as.factor(pd$Age)
pd$Pmi = as.factor(pd$Pmi)
pd$Race = as.factor(pd$Race)

#' Quantile normalization
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
write.table(p, file = paste(output_folder,suffix, "_ExpQuant", ".txt",sep = ""),sep="\t",quote=F)

#' Variance Explained
var <- VarExp(p,pd,10,FALSE)
pdf("OUTPUT_BrainSeq/Variance_Explained_BrainSeq.pdf",width=8,height=6)
plotVarExp(var,"Variance Explained")
dev.off()

#' Variance Explained by Covariates at Gene Level
cl <- makeCluster(3)
registerDoParallel(cl)
form <- ~ (1|Diagnosis) + (1|Sex) + (1|Race) + RIN + (1|Age) + (1|Pmi)
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
mod <- model.matrix(~ Diagnosis+Sex+Race+Age+Pmi+RIN, data =pd)
mod0 <- model.matrix(~Sex+Race+Age+Pmi+RIN, data = pd)

#' Create SVA object with 100 permuation 
svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100,method="two-step")
pdSv <- cbind(mod,svaobj$sv)

#' Combine the original pheno data with SVs
pdSVs<- cbind(pd,svaobj$sv)

#' Check all covariates
```{r, message = FALSE}
ggpairs(pdSVs[c(1:11)], mapping = aes(color = Diagnosis)) + 
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
BrainSeq_DGE = topTable(fitEb, coef = "DiagnosisSchizo",number=nrow(p));
write.table(BrainSeq_DGE, file = paste(output_folder,suffix, "_LIMMA_DEG", ".txt",sep = ""),sep="\t",quote=F)

#' DGE table
DT::datatable(BrainSeq_DGE, options = list(pageLength = 10))

#' Save into a .RData
pdSVs=cbind(pd,svaobj$sv)
BrainSeq_DATA=list(BrainSeq_ExpQuant=p,BrainSeq_ExpCleaned_lm=pAdj,BrainSeq_pData=pdSVs,BrainSeq_Dge_All=BrainSeq_DGE)
save(BrainSeq_DATA,file="OUTPUT_BrainSeq/BrainSeq_DATA.RData")

#' sessionInfo
sessionInfo()





