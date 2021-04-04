
multi_cox<-function(covariates,data){
  
  res.cox <- coxph(Surv(survival, death) ~ `CD4 Naive` + `CD4 Memory` + `CD8 Memory`+
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
  
  result<-list(res1,res2)
  return(result)
}
