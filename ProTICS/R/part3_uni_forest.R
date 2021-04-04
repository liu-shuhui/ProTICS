# forestplot for univariable regression
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

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(suvival,death)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data1)})

univ_results <- lapply(univ_models,function(x){
  x <- summary(x)
  p.value<-format(x$wald["pvalue"], scientific = TRUE,digits = 3)
  #wald.test<-signif(x$wald["test"], digits=2)
  #beta<-signif(x$coef[1], digits=2);#coeficient beta
  HR <-round(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower <- round(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- round(x$conf.int[,"upper .95"],2)
  HR1 <- paste0(HR, " [",HR.confint.lower, "-", HR.confint.upper, "]")
  res.cox<-c(HR, HR.confint.lower,HR.confint.upper,HR1,p.value)
  names(res.cox)<-c("HR", "lower","upper","HR [95% CI for HR]","p.value")
  return(res.cox)
  #return(exp(cbind(coef(x),confint(x))))
})

univ_res <- t(as.data.frame(univ_results, check.names = F))
res1<-rbind(c("Immune cells","HR 95% CI","P.value"),
            cbind(colnames(cell_new),format(univ_res[,c(4,5)],scientific = TRUE,digits = 3)))
res2<-data.table(rbind(c(NA,NA,NA),univ_res[,c(1,2,3)]))


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

