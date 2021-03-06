# ProTICS reveals prognostic impact of immune cell types in different molecular subtypes

ProTICS is a pipeline consisted of three parts and the three parts have their respective goals. 
The implementation of the latter parts depends on the results of the previous parts.

## Setup the environment

```{r data}
# Please install the following packages
library(data.table)
library(dplyr)
library(rTensor)
library(nnTensor)
library(survival)
library(survminer)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(grDevices)
library(pheatmap)
library(forestplot)
```

## How to use ProTICS
run Demo.R

### Give an example using a small dataset

#### Input data 
The data is in the ProTICS/Data folder.
```{r data}
# read gene expression data
data1<-fread(file = "./Data/data1.txt",header = T)
```
<div align=center><img width="550" src="https://user-images.githubusercontent.com/80741925/113571091-4850d200-9648-11eb-8fcc-eb88565797d8.png"/></div>

```{r data}
# read DNA methylation data
data2<-fread(file = "./Data/data2.txt",header = T)
```
<div align=center><img width="550" src="https://user-images.githubusercontent.com/80741925/113571112-5272d080-9648-11eb-8823-4b2ef0e70a34.png"/></div>


#### Results of part 1
Molecular subtypes discovery by running NTD method. Here, patients was clustered into two cancer subtypes.

```{r data}
## k=2 is an example. 
Subtype= NTD_subtyping(data1,data2,k=2, n=100) 
```
```{r data}
# Visualize the overall survival analysis of the two cancer subtypes
ggsurvplot(survival_out, data = survivaldata, risk.table = T,xlab="Survival time/day", ylab="Survival rate")
```
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113572452-d0d07200-964a-11eb-91e5-f07d19b9afbb.png"/></div>

#### Results of part 2

Differential expression (DE) analysis of the signature genes between the two cancer subtypes
```{r data}
GS<-subtypes_DEA(Surv,seqd)
```
Visualize the heatmaps of the selected DE genes 
```{r data}
pheatmap(data,cluster_rows=T, color = colorRampPalette(c( "#0077FF","#FFEEFF","#FF7700"))(1000),
         cluster_cols=F,show_rownames = TRUE,show_colnames=F, annotation=anno_c,annotation_legend=TRUE,main="dataset")
```
![image](https://user-images.githubusercontent.com/80741925/113573101-09bd1680-964c-11eb-989e-fc04d041c47b.png)

#### Results of part 3
```{r data}
# Distribution of proportion for the 10 immune cell types in different molecular subtypes
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
```
![image](https://user-images.githubusercontent.com/80741925/113576769-bbf7dc80-9652-11eb-95ca-73ad935711d2.png)

```{r data}
# Analysis of the prognosis of single immune cell type in subtypes 1 by using univariate cox regression
source("./R/functions/uni_cox.R")
result<-uni_cox(covariates,data1)
```
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113577769-599fdb80-9654-11eb-941a-8a54ab198fe2.png"/></div>

```{r data}
# Analysis of the prognosis of the 10 immune cell type in subtypes 1 by using multivariate cox regression
source("./R/functions/multi_cox.R")
result<-multi_cox(covariates,data1)
```
<div align=center><img width="500" src="https://user-images.githubusercontent.com/80741925/113577785-615f8000-9654-11eb-937b-4547429b8af7.png"/></div>

## How to Cite ProTICS

Shuhui Liu, Yupei Zhang, Xuequn Shang* and Zhanglei Zhang*. ProTICS reveals prognostic impact of immune cell types in different molecular subtypes, 2021, Briefing in Bioinformatics.

If you have questions on our codes, please contact me at: lsh@mail.nwpu.edu.cn
