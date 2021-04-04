##  Differential expression Analysis
subtypes_DEA <- function(Surv,seqd){
  ### define a 3-mode tensor

  group<-factor(Surv$Subtype,c("1","2"),c("Subtype_1","Subtype_2"))
  design<-model.matrix(~0+group)
  colnames(design)<-c("Subtype_1","Subtype_2")
  #y <- cpm(seqd[,-1],log = TRUE)
  y <- voom(seqd[,-1], design, plot = F)
  fit <- lmFit(y, design)
  contr <- makeContrasts(Subtype_1-Subtype_2, levels = design)
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  res <- topTable(tmp, sort.by = "P", n = Inf)
  rownames(res)<-seqd$symbol[as.numeric(rownames(res))]

  T<-res[which(abs(res$logFC)>=1 & (res$adj.P.Val < 1e-2)),]
  if (dim(T)[1]<=20)
  {GS<-rownames(T)} else {
    T<-cbind(rownames(T),T)
    colnames(T)[1]<-c("Genes")
    res<-arrange(T,desc(abs(T$logFC)))
    GS<-as.character(res[1:20,1])
  }

  return(GS)
}


