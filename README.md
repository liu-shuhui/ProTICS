# ProTICS reveals prognostic impact of immune cell types in different molecular subtypes

ProTICS is a pipeline consisted of three parts and the three parts have their respective goals. 
The implementation of the latter parts depends on the results of the previous parts.

# How to use ProTICS
run demo.R

# Input data 
We give an example using a small dataset including 200 ptients and 2000 genes. 

#Gene expression

data1<-fread(file = "./Data/data1.txt",header = T)
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113571091-4850d200-9648-11eb-8fcc-eb88565797d8.png"/></div>

#DNA methylation

data2<-fread(file = "./Data/data2.txt",header = T)
<div align=centre><img width="500" src="https://user-images.githubusercontent.com/80741925/113571112-5272d080-9648-11eb-8823-4b2ef0e70a34.png"/></div>

# part 1
#run NTD method for clustering patients into two cancer subtypes.

Subtype= NTD_subtyping(data1,data2,k=2, n=100)

#Visualize the overall survival analysis of the two cancer subtypes

ggsurvplot(survival_out, data = survivaldata, risk.table = T,xlab="Survival time/day", ylab="Survival rate")
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113572452-d0d07200-964a-11eb-91e5-f07d19b9afbb.png"/></div>

# part 2

# Part 3


