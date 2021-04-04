
###  comparison with other data-integration method for discovering cancer subtypes
library(data.table)
library(CancerSubtypes)
library("RTCGA.mRNA")
library(survival)

# load data
GE<-fread(file = "gene_expression.txt",header = T)
DM<-fread(file = "DNA_methylation.txt",header = T)
DAT=list(GeneExp=GE[,-1],GeneMehty=DM[,-1])
survdata<-fread(file = "clinical_data.txt",header = T) #include patient codes, gender, survival status and overall survival time

# iclustering method

data1=FSbyVar(GE[,-1], cut.type="topk",value=80)
data2=FSbyVar(DM[,-1], cut.type="topk",value=50)
DAT=list(GeneExp=data1,Methy=data2)
result=ExecuteiCluster(datasets=DAT, k=2, lambda=list(0.44,0.33,0.28))

survdata1<-cbind(survdata,result$group)
colnames(survdata1)<-c("patients","death","survival","Labels")
survdiff(Surv(survival,death)~Labels, data=survdata1)



# SNF
result=ExecuteSNF(DAT, clusterNum=2, K=20, alpha=0.5)
#sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)

survdata1<-cbind(survdata,result$group)
colnames(survdata1)<-c("patients","death","survival","Labels")
survdiff(Surv(survival,death)~Labels, data=survdata1)


#######CC
result=ExecuteCC(clusterNum=2,d=DAT,maxK=10,clusterAlg="hc",distance="pearson")

survdata1<-cbind(survdata,result$group)
colnames(survdata1)<-c("patients","death","survival","Labels")
survdiff(Surv(survival,death)~Labels, data=survdata1)


