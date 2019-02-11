# Cross validation for NeuN
# Slow processing
source("Utility_Functions.R")
library(ggplot2)
library(ggpubr)
load("OUTPUT_NeuN/NEUN_EXON_DATA.RData")

# Load the data and the demographic data
exp <- NEUN_EXON_DATA$NeuN_ExpQuant
pheno <- NEUN_EXON_DATA$NeuN_pData

# Run limma subsampling for 20 subjects and with 200 bootstrap
# Sample = n of sample to subset
# B = number of bootstrap
looDGE <- limmaLOO(exp,pheno,Sample=20,B=200) 
save(looDGE,file = "OUTPUT_NeuN/looDGE_NeuN.RData")

# Run limma permuting the expression 200 times
# nPERM = number of permutation
permDGE <- limmaPERM(exp,pheno,nPERM=200)
save(permDGE,file = "OUTPUT_NeuN/permDGE_NeuN.RData")

# Compare with observed estimate (logFC) with LOO estimates
obs <- NEUN_EXON_DATA$NeuN_Dge_All
obs <- obs[order(rownames(obs)),]

# Extract Estimate from loo method
extracted_cols <- lapply(looDGE, function(x) x[1])
tmp <- do.call("cbind", extracted_cols) 

# Scatter of LOO estimate and observed
newdata <- data.frame(Obs = obs$logFC, LOO = rowMeans(tmp))
pdf("OUTPUT_NeuN/LOO_Analysis_NeuN.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "LOO",
   color = "cyan4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
xlab("log2(FC) - Observed")+ 
ylab("log2(FC) - LOO")+
theme(legend.position="none")
dev.off()

# Doing it for only the DGE
dge <- obs[obs$adj.P.Val < 0.01,]

# Extract Estimate from loo method
extracted_cols <- lapply(looDGE, function(x) x[1])
tmp <- do.call("cbind", extracted_cols) 
tmp <- tmp[rownames(tmp)%in%rownames(dge),]

# Scatter of LOO estimate and observed
newdata <- data.frame(Obs = dge$logFC, LOO = rowMeans(tmp))
pdf("OUTPUT_NeuN/LOO_Analysis_NeuN_FDR001.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "LOO",
   color = "cyan4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
xlab("log2(FC) - Observed")+ 
ylab("log2(FC) - LOO")+
theme(legend.position="none")
dev.off()



# Compare with observed estimate (logFC) with Permuted estimates
extracted_cols <- lapply(permDGE, function(x) x[1])
tmp <- do.call("cbind", extracted_cols) 

# Scatter of LOO estimate and observed

newdata <- data.frame(Obs = obs$logFC, PERM = rowMeans(tmp))
pdf("OUTPUT_NeuN/PERM_Analysis_NeuN.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "PERM",
   color = "cyan4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
xlab("log2(FC) - Observed")+ 
ylab("log2(FC) - Permutation")+
theme(legend.position="none")
dev.off()

# Only significant 
dge <- obs[obs$adj.P.Val < 0.01,]

# Extract Estimate from loo method
extracted_cols <- lapply(permDGE, function(x) x[1])
tmp <- do.call("cbind", extracted_cols) 
tmp <- tmp[rownames(tmp)%in%rownames(dge),]

# Scatter of LOO estimate and observed
newdata <- data.frame(Obs = dge$logFC, PERM = rowMeans(tmp))
pdf("OUTPUT_NeuN/PERM_Analysis_NeuN_FDR001.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "PERM",
   color = "cyan4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
xlab("log2(FC) - Observed")+ 
ylab("log2(FC) - Permutation")+
theme(legend.position="none")
dev.off()

