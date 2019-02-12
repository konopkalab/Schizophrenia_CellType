suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(UpSetR))

file_names=as.list(dir(pattern=".txt"))
l <- lapply(file_names,read.table,sep="\t")

brainseq <- as.data.frame(l[[1]])
cmc <- as.data.frame(l[[2]])
neun <- as.data.frame(l[[3]])
olig <- as.data.frame(l[[4]])

cmc$Rows <- rownames(cmc)
brainseq$Rows <- rownames(brainseq)
neun$Rows <- rownames(neun)

files <- list(cmc,brainseq,neun)
NeuN_DB <- Reduce(function(x, y) {
     merge(x, y, all=FALSE, by="Rows")
}, files)

colnames(NeuN_DB) <- gsub("\\.x","_CMC",colnames(NeuN_DB))
colnames(NeuN_DB) <- gsub("\\.y","_BrainSeq",colnames(NeuN_DB))
colnames(NeuN_DB)[14:19] <- paste(colnames(NeuN_DB)[14:19],"_NeuN",sep="")

write.table(NeuN_DB,"NeuN_CMC_BrainSeq_DB.txt",sep="\t",quote=F)
write.xlsx(NeuN_DB, file="Comparative_DB.xlsx",sheetName = "NeuN",row.names=FALSE, showNA=FALSE)

# Add Olig
olig$Rows <- rownames(olig)

files <- list(cmc,brainseq,olig)
OLIG2_DB <- Reduce(function(x, y) {
     merge(x, y, all=FALSE, by="Rows")
}, files)

colnames(OLIG2_DB) <- gsub("\\.x","_CMC",colnames(OLIG2_DB))
colnames(OLIG2_DB) <- gsub("\\.y","_BrainSeq",colnames(OLIG2_DB))
colnames(OLIG2_DB)[14:19] <- paste(colnames(OLIG2_DB)[14:19],"_OLIG2",sep="")

write.table(OLIG2_DB,"OLIG2_CMC_BrainSeq_DB.txt",sep="\t",quote=F)
write.xlsx(OLIG2_DB, file="Comparative_DB.xlsx",sheetName = "OLIG2",row.names=FALSE, showNA=FALSE,append=TRUE)

#########################################
## Set of genes significant (unadj P)  ##
#########################################
res <- NeuN_DB[NeuN_DB$P.Value_CMC < 0.05 & NeuN_DB$P.Value_BrainSeq < 0.05 & NeuN_DB$P.Value_NeuN < 0.05, ] #Select genes significant in all the three data
write.xlsx(res, file="Comparative_DB_sign.xlsx",sheetName = "NeuN Nominal P",row.names=FALSE, showNA=FALSE,append=TRUE)

a <- ggscatter(res, x = "logFC_CMC", y = "logFC_NeuN",
   color = "black",size = 0.5,
   add = "reg.line",  # Add regressin line
   add.params = list(color = "steelblue", fill = "lightgrey"), 
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
xlab("CMC log2(Fold Change)")+ 
ylab("NeuN log2(Fold Change)")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
xlim(-0.2,+0.2)+
ylim(-0.65,+0.65)

b <- ggscatter(res, x = "logFC_BrainSeq", y = "logFC_NeuN",
   color = "black",size = 0.5,
   add = "reg.line",  # Add regressin line
   add.params = list(color = "steelblue", fill = "lightgrey"), 
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
xlab("BrainSeq log2(Fold Change)")+ 
ylab("NeuN log2(Fold Change)")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
xlim(-0.2,+0.2)+
ylim(-0.65,+0.65)

# OLIG2
res <- OLIG2_DB[OLIG2_DB$P.Value_CMC < 0.05 & OLIG2_DB$P.Value_BrainSeq < 0.05 & OLIG2_DB$P.Value_OLIG2 < 0.05, ] #Select genes significant in all the three data
write.xlsx(res, file="Comparative_DB_sign.xlsx",sheetName = "OLIG2 Nominal P",row.names=FALSE, showNA=FALSE,append=TRUE)

c <- ggscatter(res, x = "logFC_CMC", y = "logFC_OLIG2",
   color = "black",size = 0.5,
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgrey"), 
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
xlab("CMC log2(Fold Change)")+ 
ylab("OLIG2 log2(Fold Change)")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
xlim(-0.2,+0.2)+
ylim(-1.5,+1.5)


d <- ggscatter(res, x = "logFC_BrainSeq", y = "logFC_OLIG2",
   color = "black",size = 0.5,
   add = "reg.line",  # Add regressin line
   add.params = list(color = "magenta", fill = "lightgrey"), 
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman",label.sep = "\n"))+
xlab("BrainSeq log2(Fold Change)")+ 
ylab("OLIG2 log2(Fold Change)")+
geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
xlim(-0.2,+0.2)+
ylim(-1.5,+1.5)

plot <- plot_grid(a,b,c,d,labels=c("A","B","C","D"), ncol = 4,nrow=1,align = "h")
save_plot("FoldChange_Comparative_Analysis.pdf", plot, ncol = 2,base_height=2.5,base_width=5)

# UpSet
suppressPackageStartupMessages(library(UpSetR))

file_names=as.list(dir(pattern="DEG.txt"))
l <- lapply(file_names,read.table,sep="\t")

brainseq <- as.data.frame(l[[1]][l[[1]]$P.Value < 0.05,])
cmc <- as.data.frame(l[[2]][l[[2]]$P.Value < 0.05,])
neun <- as.data.frame(l[[3]][l[[3]]$P.Value < 0.05,])
olig <- as.data.frame(l[[4]][l[[4]]$P.Value < 0.05,])

brainseq.tmp <- data.frame(Gene = rownames(brainseq),Class = rep("BrainSeq",nrow(brainseq)))
cmc.tmp <- data.frame(Gene = rownames(cmc),Class = rep("CMC",nrow(cmc)))
neun.tmp <- data.frame(Gene = rownames(neun),Class = rep("NeuN",nrow(neun)))
olig.tmp <- data.frame(Gene = rownames(olig),Class = rep("OLIG2",nrow(olig)))
df <- rbind(brainseq.tmp,cmc.tmp,neun.tmp,olig.tmp)

final.l <- split(as.character(df$Gene),df$Class)
Class <- names(final.l)
ToTGene <- as.numeric(sapply(final.l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("UpSet_NominalP.pdf",width=6,height=5,useDingbats=FALSE)
upset(fromList(final.l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), list(type = "matrix_rows", 
    column = "sets", colors = c(NeuN = "cyan4", OLIG2 = "magenta4", CMC = "grey60",BrainSeq="steelblue"), 
    alpha = 0.5))))
dev.off()
