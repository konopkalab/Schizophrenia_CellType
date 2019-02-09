# Cross validation for NeuN
# Slow processing
source("Utility_Functions.R")
load("OUTPUT_OLIG2/OLIG2_EXON_DATA.RData")

# Load the data and the demographic data
exp <- OLIG2_EXON_DATA$OLIG2_ExpQuant
pheno <- OLIG2_EXON_DATA$OLIG2_pData

# Run limma subsampling for 15 subjects and with 200 bootstrap
# Sample = n of sample to subset
# B = number of bootstrap
looDGE <- limmaLOO(exp,pheno,Sample=15,B=200) 
save(looDGE,file = "OUTPUT_OLIG2/looDGE_OLIG2.RData")

# Run limma permuting the expression 200 times
# nPERM = number of permutation
permDGE <- limmaPERM(exp,pheno,nPERM=200)
save(permDGE,file = "OUTPUT_OLIG2/permDGE_OLIG2.RData")

# Compare with observed estimate (logFC) with LOO estimates
obs <- OLIG2_EXON_DATA$OLIG2_Dge_All
obs <- obs[order(rownames(obs)),]

# Extract Estimate from loo method
extracted_cols <- lapply(looDGE, function(x) x[1])
tmp <- do.call("cbind", extracted_cols) 

# Scatter of LOO estimate and observed
newdata <- data.frame(Obs = obs$logFC, LOO = rowMeans(tmp))
pdf("OUTPUT_OLIG2/LOO_Analysis_OLIG2.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "LOO",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
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
pdf("OUTPUT_OLIG2/LOO_Analysis_OLIG2_FDR001.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "LOO",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
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
pdf("OUTPUT_OLIG2/PERM_Analysis_OLIG2.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "PERM",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
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
pdf("OUTPUT_OLIG2/PERM_Analysis_OLIG2_FDR001.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(newdata, x = "Obs", y = "PERM",
   color = "magenta4", shape = 21, size = 1, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
xlab("log2(FC) - Observed")+ 
ylab("log2(FC) - Permutation")+
theme(legend.position="none")
dev.off()

