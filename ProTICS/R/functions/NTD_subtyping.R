# NTD_subtyping
#
# This function is used to integrate multi-omics and
# perform Non-negative Tucker Decomposition Algorithms
# then assign patients into different groups by matrice_B.

NTD_subtyping <- function(data1,data2,k,n){
  ### define a 3-mode tensor
  
  arr <- array(0,dim = c(dim(data1[,-1]),2)) # rows: genes; columns: patiens
  arrT <- as.tensor(arr)
  
  arrT[,,1] <- unlist(normalization(data1[,-1]))
  arrT[,,2] <- unlist(normalization(data2[,-1]))

  ##k:the number of subtypes; n:The number of interation step (Default: 100)
  
  output <- NTD(arrT, rank=c(k, k, k),num.iter=n)  # NTD function from nnTensor

  matrice_B<-t(output$A[[2]]) # matrice_B saved the latent factors information of patiens
  group<-max.col(matrice_B) # subtypes information
  return(group)
}
