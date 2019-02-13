###############################################################################################
#Steps to run DSS to identify differentially methylated CpGs in OLIG2 vs NeuN among control individuals (25 NeuN and 20 OLIG2 samples)
###############################################################################################
library(bsseq)
library(DSS)


#read bsseq object with all sindividual and all CpGs for controls
data1<-readRDS(paste0("bsseq.Rds")) 

#read txt file with covariates (see values in Supplementary Tables). The file contains as many rows as samples and 13 columns:
data2<-read.table("covariates_WGBS.txt",header=T)
#col1:Sample_ID
#col2:Sample_ID_2 (not necessary in this analysis)
#col3:Diagnosis (not necessary in this analysis)
#col4:Cell_type 
#col5:Sex
#col6:Age (not necessary in this analysis)
#col7:AgeClass
#col8:BrainBank
#col9:Pmi (not necessary in this analysis)
#col10:PmiClass
#col11:Hemisphere
#col12:Conversion_rates
#col13: Sample_ID_PCs

#read txt file with 10 genetic PCs. The file contains as rows as individuals (54) and 10 columns (each PC)
data3<-read.table("top_10PC_scores.txt",header=T)

#1-put covariate file in the same order as the individuals in the bsseq object
match1<-match(colnames(data1),data2[,1])
data4n<-data2[match1,]

#2-one sample might have NA conversion rate, put mean instead
if(sum(is.na(data4n[,12]))>=1) {   data4n[is.na(data4n[,12]),12]<-mean(na.omit(data4n[,12])) } 

#3-put PCs for each individual in same order as covariates
match2<-match(data4n[,13],rownames(data3))
data3n<-data3[match2,]


#4-run DSS for diagnosis using covariates + 10 genetic PCs in the model
design<-data.frame(as.factor(data4n[,4]),as.factor(data4n[,5]),as.factor(data4n[,7]),as.factor(data4n[,8])
                  ,as.factor(data4n[,10]),as.factor(data4n[,11]),as.numeric(data4n[,12]),data3n)

colnames(design)<-c("Cell_type","sex","age","bank","pmi","hemi","conv",paste0("pc",1:10))


DMLfit = DMLfit.multiFactor(data1, design, ~Cell_type+sex+age+bank+pmi+hemi+conv+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10)
res = DMLtest.multiFactor(DMLfit, coef=2)
save.image(paste("DSS_celltype_10geneticPCs.Rds",sep=''))
write.csv(res,paste0("DSS_celltype_10geneticPCs.csv"),quote=F,row.names=F)


###############################################################################################
#now identify DMRs using DSS with DMPs at bonferroni P<0.05 threshold
####################################################################################
total <- dim(res)[1] - length(which(is.na(res[,4])))
res.dmr <- callDMR(res,p.threshold=(0.05/total),minCG = 5,dis.merge = 100)
write.table(res.dmr, "DSS_celltype_DMRs.txt")
###################################################################################


###############################################################################################
#Run LM on the same model
###############################################################################################
meth.df <- getMeth(data1,type="raw")
lm.res <- matrix(nrow=dim(meth.df)[1],ncol=4)
for(i in 1:dim(meth.df)[1]) {
	try(cur.res <- lm(as.numeric(meth.df[i,])~design$Cell_type+design$sex+design$age+design$bank+design$pmi+design$hemi+design$conv+design$pc1+design$pc2+design$pc3+design$pc4+design$pc5+design$pc6+design$pc7+design$pc8+design$pc9+design$pc10))
    try(lm.res[i,1]<- as.numeric(summary(cur.res)$coefficients[2,1]))
    try(lm.res[i,2]<- as.numeric(summary(cur.res)$coefficients[2,2]))
    try(lm.res[i,3]<- as.numeric(summary(cur.res)$coefficients[2,3]))
    try(lm.res[i,4]<- as.numeric(summary(cur.res)$coefficients[2,4]))
    if(i%%20==0){print(i)}
}
colnames(lm.res) <- c("estimate","se","tvalue","pvalue")
write.csv(lm.res,paste0("DSS_celltype_10geneticPCs.csv"),quote=F,row.names=F)


