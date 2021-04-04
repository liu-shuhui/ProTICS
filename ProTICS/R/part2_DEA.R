## this Differential expression Analysis and Plot heatmaps for the selected top genes
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(grDevices)
library(data.table)
library(dplyr)

## Read the counts
sig_expr <- fread("gene_exprssion.txt",sep = "\t",header = TRUE) ###rows are the signature genes
survival_data <- fread("survival_data.txt", sep = "\t",header = TRUE) ### subtype_label can be saved from part1.R
subtypes<-survival_data$Labels

####################################################S1 VS S2 ####################################
ID1<- which(subtypes==1 | subtypes==2)
Surv1<-survival_data[ID1,]
seqd1<-select(sig_expr,c(colnames(sig_expr)[1],Surv1$patient_id))

group1<-factor(Surv1$Labels,c("1","2"),c("Subtype_1","Subtype_2"))
design1<-model.matrix(~0+group1)
colnames(design1)<-c("Subtype_1","Subtype_2")
#y <- cpm(seqd1[,-1],log = TRUE)
y <- voom(seqd1[,-1], design1, plot = F)
fit1 <- lmFit(y, design1)
contr1 <- makeContrasts(Subtype_1-Subtype_2, levels = design1)
tmp1 <- contrasts.fit(fit1, contr1)
tmp1 <- eBayes(tmp1)
res1 <- topTable(tmp1, sort.by = "P", n = Inf)
rownames(res1)<-seqd1$symbol[as.numeric(rownames(res1))]

T1<-res1[which(abs(res1$logFC)>=1 & (res1$adj.P.Val < 1e-2)),]
if (dim(T1)[1]<=20)
{GS1<-rownames(T1)} else {
  T1<-cbind(rownames(T1),T1)
  colnames(T1)[1]<-c("Genes")
  res1<-arrange(T1,desc(abs(T1$logFC)))
  GS1<-as.character(res1[1:20,1])
}
UCEC<-GS1
# GS1<-rownames(T1)


###################################################heatmaps #####################################

sig_expr<-sig_expr[is.element(sig_expr$symbol,GS1),]

IDD<-c(which(subtypes==1),which(subtypes==2))
survd_new<-survival_data[IDD,]

sigdata<-select(sig_expr,c(colnames(sig_expr)[1],survd_new$patient_id))

anno_c<-data.frame(Types = factor(survd_new$Labels,c("1","2"),c("S1","S2")))
colnames(anno_c)<-c("  ")
row.names(anno_c)<-survd_new$patient_id

data<-normalization(sigdata[,-1])

library(pheatmap)
rownames(data)<-sigdata$symbol
pheatmap(data,cluster_rows=T,
         color = colorRampPalette(c( "#0077FF","#FFEEFF","#FF7700"))(1000),
         cluster_cols=F,show_rownames = TRUE,show_colnames=F,
         annotation=anno_c,annotation_legend=TRUE,main="dataset")
############################################################# Heatmaps

