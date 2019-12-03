## Creating a basic PCA plot from a vcf file 
This document describes an approach to visualise the population genetic structure using a Principle Components Analysis (PCA). We begin by loading the relaveant libraries, setting the working directory, then reading in vcf files and convert to gds.   

```
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(vcfR)
library(hierfstat)
library(adegenet)
```
Load the vcf files and convert to ```gds```using the ```SNPRelate::snpgdsVCF2GDS```
```
file1 = "file1.vcf"
file2 = "file2.vcf"
file3 = "file3.vcf"

snpgdsVCF2GDS(file1, "original.gds")
snpgdsVCF2GDS(file2, "LDpruned.gds")
snpgdsVCF2GDS(file3, "filtered.gds")

o.genofile <- snpgdsOpen("original.gds")
p.genofile <- snpgdsOpen("LDpruned.gds")
f.genofile <- snpgdsOpen("filtered.gds")

```

Then we perform a PCA on the datasets. We make an object (i.e. o.pc.percent) which is the proportion of variance explained by the princple components. We use this later to label our PCA axis. We also create a dataframe with the first two PCs.
```
o.pca<-snpgdsPCA(o.genofile,autosome.only=FALSE)
p.pca<-snpgdsPCA(p.genofile,autosome.only=FALSE)
f.pca<-snpgdsPCA(f.genofile,autosome.only=FALSE)

o.pc.percent <- o.pca$varprop*100
p.pc.percent <- p.pca$varprop*100
f.pc.percent <- f.pca$varprop*100

o.Xlab<-paste(paste("PCA1 (", round(o.pc.percent, 2)[1], sep=""),"%)", sep="")
o.Ylab<-paste(paste("PCA2 (", round(o.pc.percent, 2)[2], sep=""),"%)", sep="")

p.Xlab<-paste(paste("PCA1 (", round(p.pc.percent, 2)[1], sep=""),"%)", sep="")
p.Ylab<-paste(paste("PCA2 (", round(p.pc.percent, 2)[2], sep=""),"%)", sep="")

f.Xlab<-paste(paste("PCA1 (", round(f.pc.percent, 2)[1], sep=""),"%)", sep="")
f.Ylab<-paste(paste("PCA2 (", round(f.pc.percent, 2)[2], sep=""),"%)", sep="")

data.o <-data.frame(pca1=o.pca$eigenvect[,1],pca2=o.pca$eigenvect[,2], sample=o.pca$sample.id)
data.p <-data.frame(pca1=p.pca$eigenvect[,1],pca2=p.pca$eigenvect[,2], sample=p.pca$sample.id)
data.f <-data.frame(pca1=f.pca$eigenvect[,1],pca2=f.pca$eigenvect[,2], sample=f.pca$sample.id)
```
## Adding metadata - phenotypic or population based information

We will add metadata to the data frame we created. Load the relevant txt file that came with the data.
If this file does not work  - can you figure out why and can you figure out how to read this into R. Try solving it with R or BASH, not excel. 
Once you have solved this problem - we will them merge the metadata file with our data frame containing the PCs. To merge we need two data frames with a matching column - some or all of the values in the matching column must be identical in both data frames and the column must have the same header - in this case we are matching individual sample ID between the PC data frame and the metadata file. 
You can check the names of the dataframes using the ```names()```function. You may need to change the column name in the metadata dataframe. Then use merge as follows

```
o.db <- droplevels(merge(data.o, metadata))
p.db <- droplevels(merge(data.p, metadata))
f.db <- droplevels(merge(data.f, metadata))
```
How many rows and columns (dimensions) does the new dataframe have? How does this compare with the dimensions of data.o?
The all.x=TRUE tells R to add rows to the new dataframe with NAs is if there are values in x that are not in y. Do we need to do this?

## Create a basic PCA plots for the three datasets

```
ggplot(o.db, aes(x=pca1, y=pca2))+
  geom_point()+
  xlab(o.Xlab)+
  ylab(o.Ylab)

ggplot(p.db, aes(x=pca1, y=pca2))+
  geom_point()+
  xlab(p.Xlab)+
  ylab(p.Ylab)

ggplot(f.db, aes(x=pca1, y=pca2))+
  geom_point()+
  xlab(f.Xlab)+
  ylab(f.Ylab)

```
### Going further
We could add colours to the plots, e.g. if your data frame has population information you can colour the dots accordingly. NB if there are many levels of the grouping factor then this is not advisable. To see how many levels there are to a factor use ```length(levels(dataframe$colname))```

Example of how to add colour to the plot
```
ggplot(f.db, aes(x=pca1, y=pca2, colour=Population))+
  geom_point()+
  xlab(f.Xlab)+
  ylab(f.Ylab)
```
The website below is a tutorial for PCAplots by Boqiang Hu and is good place to start.
http://huboqiang.cn/2016/03/03/RscatterPlotPCA
