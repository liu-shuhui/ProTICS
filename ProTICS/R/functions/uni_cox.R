
uni_cox<-function(covariates,data){
  
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(survival,death)~', x)))
  
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

  result<-list(res1,res2)
  return(result)
}
