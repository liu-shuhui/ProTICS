# forestplot for multivariable regression
###
### run once for each subtype
### change Labels==i for subtype i

library(forestplot)
library(data.table)
library(survival)
library("survminer")
library(tableone)
library(dplyr)

survdata<-fread("survdata.txt")
survd<-survdata[which(Labels==1),] ## This runs for subtype 1;
survd$suvival<-scale(survd$suvival,center = FALSE, scale = TRUE)

cell<-fread(file = "CIBER_User_Allsample_Result.txt", sep = "\t",header = T)
cell_n<-cell[,-c(2,5,8,14,16:18)] ## remove the clolumns which are not cell types.
ID1<-sapply(survd$patient_id,function(x) which(cell_n$Mixture==x))
cell_new<-cell_n[ID1,-1]

cutoff<-as.matrix(apply(cell_new,2,median)) #### calculate the median value for each cell type in a subtype;

A<-t(replicate(dim(cell_new)[1],cutoff[,1]))
C<-as.matrix(cell_new>A)
C[C==TRUE]<-1

data1<-cbind(survd,C)
covariates<-c("`CD4 Naive`","`CD4 Memory`","`CD8 Memory`",
              "`CD8 Effector`", "`Th cell`", "`Monocytes CD16`",
              "`Monocytes CD14`","DC","pDC","Plasma")

res.cox <- coxph(Surv(suvival, death) ~ `CD4 Naive` + `CD4 Memory` + `CD8 Memory`+
                   `CD8 Effector`+`Th cell`+`Monocytes CD16`+`Monocytes CD14`+DC+
                   pDC+Plasma, data =  data1)
#summary(res.cox)

multi_res <- summary(res.cox)
res1 <- cbind(colnames(cell_new),multi_res[["coefficients"]][,c(2,5)])
res2<-multi_res[["conf.int"]][,-2]

HR <-round(res2[,1], digits=2);#exp(beta)
HR.confint.lower <- round(res2[,2], 2)
HR.confint.upper <- round(res2[,3],2)
res1[,2] <- paste0(HR, " [",HR.confint.lower, "-", HR.confint.upper, "]")
res1[,3]<-format(as.numeric(res1[,3]), scientific = TRUE, digits = 2)
res1 <- rbind(c("Immune cells","HR 95% CI","P.value"),res1)

res2<-data.table(rbind(c(NA,NA,NA),res2))
colnames(res2)<-c("HR", "lower","upper")
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



