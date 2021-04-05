library(data.table)
library(dplyr)
library(rTensor)
library(nnTensor)
library(survival)
library(survminer)

###################### Part1 ######################
setwd("C:/Users/Yupei Zhang/Desktop/lsh/bib/minor revision/ProTICS/")

data1<-fread(file = "./Data/data1.txt",header = T)
data2<-fread(file = "./Data/data2.txt",header = T)
clinicdata<-fread(file ="./Data/clinic_Data.txt",header = T)
colnames(clinicdata)<-c("patient_id", "death", "survival")

source("./R/functions/normalization.R")
source("./R/functions/NTD_subtyping.R")

Subtype= NTD_subtyping(data1,data2,k=2, n=100)

survivaldata<-cbind(clinicdata,Subtype)

write.table(survivaldata, file = "./output/overallsurvival_subtypes.txt",
            sep = "\t", col.names = T, quote = F, row.names = F)

survdiff(Surv(survival,death)~Subtype, data=survivaldata)

survival_out<-survfit(Surv(survival,death)~Subtype, data=survivaldata)

ggsurvplot(survival_out, data = survivaldata, risk.table = T,xlab="Survival time/day", ylab="Survival rate")


###################### Part2 ######################
#### Differential expression Analysis and Plot heatmaps for the selected top genes
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(grDevices)
library(data.table)
library(dplyr)

## Read the counts
sig_expr <- fread("./Data/signature_count.txt",sep = "\t",header = TRUE) ###rows are the signature genes
survival_data <- fread("./output/overallsurvival_subtypes.txt", sep = "\t",header = TRUE) ### subtype_label can be saved from part1.R
subtypes<-survival_data$Subtype

ID<- which(subtypes==1 | subtypes==2)
Surv<-survival_data[ID,]
seqd<-select(sig_expr,c(colnames(sig_expr)[1],Surv$patient_id))

source("./R/functions/subtypes_DEA.R")
GS<-subtypes_DEA(Surv,seqd)

# Heatmaps of differentially expressed genes

sig_expr<-sig_expr[is.element(sig_expr$symbol,GS),]

IDD<-c(which(subtypes==1),which(subtypes==2))
survd_new<-survival_data[IDD,]

sigdata<-select(sig_expr,c(colnames(sig_expr)[1],survd_new$patient_id))

anno_c<-data.frame(Types = factor(survd_new$Subtype,c("1","2"),c("Sub1","Sub2")))
colnames(anno_c)<-c("  ")
row.names(anno_c)<-survd_new$patient_id

source("./R/functions/normalization.R")


data<-normalization(log2(sigdata[,-1]+1))

library(pheatmap)
rownames(data)<-sigdata$symbol
pheatmap(data,cluster_rows=T,
         color = colorRampPalette(c( "#0077FF","#FFEEFF","#FF7700"))(1000),
         cluster_cols=F,show_rownames = TRUE,show_colnames=F,
         annotation=anno_c,annotation_legend=TRUE,main="dataset")


###################### Part3 ######################
library(forestplot)
library(data.table)
library(survival)
library("survminer")
library(dplyr)

survdata <- fread("./output/overallsurvival_subtypes.txt", sep = "\t",header = TRUE)
cell<-fread(file = "./Data/CellProportion.txt", sep = "\t",header = T)

## the clolumns[16:18] which are not cell types should be removed.
cell<-cell[,-c(16:18)]

id=which(apply(cell[,-1],2,var)>1e-05)+1 ### remove the clolumns whose var is much small
cell_new<-select(cell,c(colnames(cell)[c(1,id)]))

### colnames of cell types
covariates<-c("`CD4 Naive`","`CD4 Memory`","`CD8 Memory`",
              "`CD8 Effector`", "`Th cell`", "`Monocytes CD16`",
              "`Monocytes CD14`","DC","pDC","Plasma")

##### 1. plot proportion distributions of the 10 cell types in different cancer subtypes
`Cell types` = c(rep(covariates, each=length(which(survdata$Subtype==1))),
                 rep(covariates, each=length(which(survdata$Subtype==2))))
`Patient type` = c(rep(c("Subtyp1"),each=length(which(survdata$Subtype==1))*10),
                   rep(c("Subtyp2"),each=length(which(survdata$Subtype==2))*10))

ID1<-sapply(survdata$patient_id[which(survdata$Subtype==1)],
                 function(x) which(cell_new$Mixture==x))
ID2<-sapply(survdata$patient_id[which(survdata$Subtype==2)],
                 function(x) which(cell_new$Mixture==x))

`Relative proportions of the 10 immune cell types` <-c(as.vector(as.matrix(cell_new[ID1,-1])),
                                                  as.vector(as.matrix(cell_new[ID2,-1])))

data<-data.frame(`Cell types`,`Patient type`,`Relative proportions of the 10 immune cell types`)


data$Cell.type <- factor(data$Cell.type,levels=covariates,ordered = TRUE)
ggplot(data, aes(`Cell types`, y=`Relative proportions of the 10 immune cell types`, color=`Patient type`)) +
  theme(
    panel.background = element_rect(linetype = 1, colour = "white", size = 1,fill = "lightblue"),
    axis.text.x = element_text(angle = 20, hjust = 0.6,vjust = 0.75),
    plot.title = element_text(colour = "black",face = "bold",size = 12, vjust = 1),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches")
  )+
  stat_boxplot(geom ='errorbar', width = 0.8) +
  geom_boxplot(width = 0.8)
facet_grid(.~Cell.type, scales = "free_x")

##### 2. forest plot of prognostic association of immune cell types

surv_sub<-survdata[which(survdata$Subtype==1),] ## This runs for subtype 1;
surv_sub$survival<-scale(surv_sub$survival,center = FALSE, scale = TRUE)


ID<-sapply(surv_sub$patient_id, function(x) which(cell_new$Mixture==x))
cell_new<-cell_new[ID,-1]
#cell_new<-logcell<-log2(cell_new+1)
cutoff<-as.matrix(apply(cell_new,2,median))

tem<-t(replicate(dim(cell_new)[1],cutoff[,1]))
mat_bip<-as.matrix(cell_new>tem)
mat_bip[mat_bip==TRUE]<-1

data1<-cbind(surv_sub,mat_bip)

## univariate cox regression
source("./R/functions/uni_cox.R")
result<-uni_cox(covariates,data1)
res1<-result[[1]]
res2<-result[[2]]
# forest plot
forestplot(res1, mean = res2$HR, lower = res2$lower, upper = res2$upper,
           graph.pos = 2,graphwidth = unit(18,"mm"),
           hrzl_lines = list("2" = gpar(lty=2,columns=1:4)),
           is.summary = c(TRUE,rep(FALSE,10)),
           txt_gp = fpTxtGp(ticks = gpar(cex=0.8),summary = gpar(cex=0.8),cex = 0.8),
           boxsize = 0.2,
           line.margin = unit(6,"mm"),
           lineheight = unit(6,"mm"),
           col=fpColors(box="blue",line="blue",summary="blue"),
           clip = c(0,5),
           xticks = c(0, 0.5, 1, 2,3,4,5),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.12,
           colgap = unit(2,"mm"),zero = 1,
           title = "Subtype 1")


#### multivariate cox regression
source("./R/functions/multi_cox.R")
result<-multi_cox(covariates,data1)
res1<-result[[1]]
res2<-result[[2]]
### ploting
forestplot(res1, mean = res2$HR, lower = res2$lower, upper = res2$upper,
           graph.pos = 2,graphwidth = unit(18,"mm"),
           hrzl_lines = list("2" = gpar(lty=2,columns=1:4)),
           is.summary = c(TRUE,rep(FALSE,10)),
           txt_gp = fpTxtGp(ticks = gpar(cex=0.8),summary = gpar(cex=0.8),cex = 0.8),
           boxsize = 0.2,
           line.margin = unit(6,"mm"),
           lineheight = unit(6,"mm"),
           col=fpColors(box="blue",line="blue",summary="blue"),
           clip = c(0,5),
           xticks = c(0, 0.5, 1, 2,3,4,5),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.12,
           colgap = unit(2,"mm"),zero = 1,
           title = "Subtype 1")

