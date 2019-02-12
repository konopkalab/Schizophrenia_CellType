library(qdap)
library(ggpubr)
library(cowplot)
library(ggthemes)

load("CellType_DMR.RData")
dge <- read.table("CELLTYPE_LIMMA_DGE.txt")

df <- merge(dge,dml, by.x="row.names",by.y="SYMBOL",all=F)
df <- na.omit(df)
df$Class <-  genX(df$annotation, " (", ")")
df$Class <- as.factor(df$Class)
df$CellType <- ifelse(df$logFC > 0, "OLIG2","NeuN")
df$DGE <- ifelse(df$Bonf < 0.05 & abs(df$logFC) > 0.5,"DGE","NOT_DGE")
prom <- df[grep("Promoter",df$Class),]

pdf("Promoter_CellType_DGE_Methylation.pdf",width=3,height=3,useDingbats=FALSE)
ggscatter(prom, x = "areaStat", y = "logFC",
   color = "DGE",size = 0.1,shape = 21, palette=c("black","lightgrey"), # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
theme_classic()+
geom_vline(xintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5) +
geom_hline(yintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5)+
theme(legend.position=c(0.8,0.9))+
xlim(-3000,+3000) +
xlab("Promoter DMR")+ 
ylab("log2(Fold-Change)")+
ggtitle("OLIG2/NeuN")
#+ geom_text_repel(data = top_labelled2, mapping = aes(label = Row.names), size = 2,color = 'black',box.padding = unit(0.1, "lines"),point.padding = unit(0.1, "lines"))
dev.off()

pdf("CellType_DGE_Methylation.pdf",width=8,height=7,useDingbats=FALSE)
ggscatter(df, x = "areaStat", y = "logFC",
   color = "DGE",size = 0.5,shape = 21, palette=c("black","grey60"), # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "spearman", label.sep = "\n")
   )+
theme_classic()+
geom_vline(xintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5) +
geom_hline(yintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5)+
theme(legend.position="none")+
xlim(-3000,+3000) +
xlab("Promoter DMR")+ 
ylab("log2(Fold-Change)")+
ggtitle("OLIG2/NeuN")+
facet_wrap(~Class)
dev.off()