# Enrichment analysis
library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
enableWGCNAThreads()
# Load data
load("NEUN_EXON_DATA.RData")

dge <- NEUN_EXON_DATA$NeuN_Dge_All
fdr001 <- dge[dge$adj.P.Val < 0.01,]
fdr001$Class <- ifelse(fdr001$logFC > 0,"Up","Down")
fdr001$Gene <- rownames(fdr001)
rownames(fdr001) <- NULL
tmp2 <- fdr001[c(8,7)]
write.table(tmp2,"NeuN_DGE_Data.txt",quote=F,sep="\t")

# Sord DSS stats and add columns
load("NeuN_DSS_20x_08ind_Filtered.RData")
classes <- split(df.n,df.n$Class)

for (i in 1:length(classes)){
	classes[[i]] <- classes[[i]][order(classes[[i]]$dss.pvalue),]
	classes[[i]]$DML_Top1k <- c(rep("DML_Top1k",1000),rep("NOT_DML",nrow(classes[[i]])-1000))
}

# Geneset DML_Top1k
names <- names(classes)
names <- gsub("' | ","_",names)
tmp <- list()
for (i in 1:length(classes)){
	tmp[[i]] <- classes[[i]][grep("DML_Top1k",classes[[i]]$DML_Top1k),]
	tmp[[i]] <- data.frame(Gene = tmp[[i]]$SYMBOL, Class = rep(names[[i]],nrow(tmp[[i]])))
}

GeneSets_Top1k <- tmp
names(GeneSets_Top1k) <- names
save(GeneSets_Top1k,file = "GeneSets_Top1k.RData")

## Enrichment analysis
table_input = list.files(pattern = 'txt')
tab=read.table(table_input,sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION")
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("GeneSets_Top1k.RData")
#for(i in 1:length(GeneSets))
#{
# GeneSets[[i]] <- GeneSets[[i]][GeneSets[[i]]$Gene %in% tab$Gene,]
#}
GeneSets <- GeneSets_Top1k
ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]])
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Pval Adjusted
library(magrittr)
FisherAdj <- FisherP %>% 
			as.matrix %>% 
			as.vector %>% 
			p.adjust(method='BH') %>% 
			matrix(ncol=ncol(FisherP))

rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj>0.05]=1
FisherOR[FisherOR < 1]=0
pdf("Heatmap_GeneSets_ProtCod_Top1k_FDR001.pdf",width=5,height=3)
df=-log10(FisherAdj)
LabelMat = paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df, 
	xLabels = colnames(df), 
	yLabels = rownames(df), 
	colorLabels =FALSE,
	colors=colorRampPalette(c("white", "red"))(50),
	textMatrix=LabelMat, 
	setStdMargins = FALSE, 
	cex.text = 0.8)
dev.off()