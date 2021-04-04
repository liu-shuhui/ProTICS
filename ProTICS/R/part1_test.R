library(data.table)
library(dplyr)
library(rTensor)
library(nnTensor)
library(survival)
library(survminer)

########load data and normalization

setwd("set path")
load("datasets.RData")
load("clinicaldata.RData") # rows: patients; columns:("patient_id", "death", "suvival");

group= NTD_subtyping(data1,data2)
survivaldata<-cbind(clinicaldata,group)

# survival analysis for the identified cancer subtypes

#survdiff(Surv(suvival,death)~group, data=survivaldata)
survival_out<-survfit(Surv(suvival,death)~group, data=survivaldata)

# plot survival curves
ggsurvplot(survival_out, data = survivaldata, risk.table = T,xlab="Survival time/day", ylab="Survival rate")


