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
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113571112-5272d080-9648-11eb-8823-4b2ef0e70a34.png"/></div>

# part 1
#run NTD method for clustering patients into two cancer subtypes.

Subtype= NTD_subtyping(data1,data2,k=2, n=100)

#Visualize the overall survival analysis of the two cancer subtypes

ggsurvplot(survival_out, data = survivaldata, risk.table = T,xlab="Survival time/day", ylab="Survival rate")
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113572452-d0d07200-964a-11eb-91e5-f07d19b9afbb.png"/></div>

# part 2

#Differential expression (DE) analysis of the signature genes between cancer subtypes

GS<-subtypes_DEA(Surv,seqd)

#Visualize the heatmaps of the selected DE genes 

pheatmap(data,cluster_rows=T, color = colorRampPalette(c( "#0077FF","#FFEEFF","#FF7700"))(1000),
         cluster_cols=F,show_rownames = TRUE,show_colnames=F, annotation=anno_c,annotation_legend=TRUE,main="dataset")

![image](https://user-images.githubusercontent.com/80741925/113573101-09bd1680-964c-11eb-989e-fc04d041c47b.png)

# Part 3

#Distribution of proportion for the 10 immune cell types in different ccancer sybtypes

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

![image](https://user-images.githubusercontent.com/80741925/113576769-bbf7dc80-9652-11eb-95ca-73ad935711d2.png)

#Analysis of the prognosis of single immune cell type in subtypes 1 by using univariate cox regression
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113577769-599fdb80-9654-11eb-941a-8a54ab198fe2.png"/></div>

#Analysis of the prognosis of the 10 immune cell type in subtypes 1 by using multivariate cox regression
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113577785-615f8000-9654-11eb-937b-4547429b8af7.png"/></div>

