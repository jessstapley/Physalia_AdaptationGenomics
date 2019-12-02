# Investigating Population Genetic Struture in R
This exercise we invesitgte population genetic strcuture using two approaches - a distance based method (PCA) and a model based method (entropy criterion).

First, create a new directory called PopGen. In that directory create one subdirectory (data). Copy the data file to the data directory.

Open R studio, create a new project using this existing directory that you just created (PopGen). In R open a new R Markdown file. Have a look at the default file format. You can add chuncks of code to the sections between "```"  "```". You can replace the code in all sections EXCEPT this bit - DO NOT EDIT this bit.
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Get familiar with this sytax and organisation of the markdown documents. Save this document (PopGenStructure) to the directory PopGen and now you can edit this document by copying the information from this github tutporial.

Now follow this tutorial on github - copy the code from here into the markdown document in a similar fashion. To run a chunk of code in the markdwon document press the small green arrow.

# PCA analysis from vcf file
Now the R project is set up lets begin.
Read in your libraries.

```
library(LEA)
library(RColorBrewer)

```
You will need to change the name of your input.path  and out.path in the code below. First we convert the vcf to gds format for this PCA analysis
```
input.path <-  "data/FILENAME.vcf"
out.path <- "data/FILENAME.geno"

vcf2geno(input.file, output.file =out.path, force = TRUE)
o.genofile <- snpgdsOpen(paste0(out.path,".gds"))
```
Then we perform a PCA on the datasets. We make an object (i.e. o.pc.percent) which is the proportion of variance explained by the princple components. We use this later to label our PCA axis. We also create a dataframe with the first two PCs.

```
o.pca<-snpgdsPCA(o.genofile,autosome.only=FALSE)
o.pc.percent <- o.pca$varprop*100
o.Xlab<-paste(paste("PCA1 (", round(o.pc.percent, 2)[1], sep=""),"%)", sep="")
o.Ylab<-paste(paste("PCA2 (", round(o.pc.percent, 2)[2], sep=""),"%)", sep="")

data.o <-data.frame(pca1=o.pca$eigenvect[,1],pca2=o.pca$eigenvect[,2], sample=o.pca$sample.id)

ggplot(o.db, aes(x=pca1, y=pca2))+
  geom_point()+
  xlab(o.Xlab)+
  ylab(o.Ylab)
  
```
Looking at the number of principle components that explain variation in the genomic data we can get an idea about the likely number of ancestral populations. This  can help us choose a range of K for the next anaysis step.

Here we perform a pca and then we use the Tracy-Widom test for each eigenvalue to investigate the number of genetic clusters in the data. We display the p-values for the Tracy-Widom test for 1-10 eigenvalues. The 'knee' in the plot indicates the number of significant components (K) in the data. The number of genetic clusters is K+1.

```
pca1R <- pca(out.path, scale =TRUE)
tw <- tracy.widom(pca1R)
tw$pvalues[1:10]
plot(tw$percentage, pch=19, typ="b") 

```
Q1. How many populations do you think there are in your data?

# Looking at population structure 

Next we will investigate the population structure more throughly. The ```snmf``` function estimates ancestry coefficients similar to the commonly used programs (STRUCTURE  and ADMIXTURE). The function estimates entropy criterion that evaluates the quality of fit of the statistical model to the data using a cross-validation technique and can help you choose the number of population clusters there are in the data. We are performing ten replicates``` rep=10```. This fuction can take some time but in comparison to running STRUCTURE it is much faster. We explore a range of values of Ks that span the true number of ancestral populations. You will need to choose a range of Ks based on your answer to Q1. E.g. if the K =3, then you could choose K=1:6.

```
obj.snmf <- snmf("file.lfmm", K=1:6, rep=10, entropy=T, ploidy =2, project ="new", iterations=10000)
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
summary(obj.snmf)
```
The plot shows the cross-entropy criterion plotted against different values of K. The lowest value of cross-entropy criterion indicates the most likely number of ancestral populations in this data. Determining the correct value of K can be difficult under some demogrpahies. If the line wisas declining this would indicate finer popuation structure possibly due to isolation by distance. The summary option shows you the min, mena and max of the cross entropy values for each K. 

Q2. How many ancetral populations are evident in your data?

We will now plot the ancestry coefficients for each individual in a bar plot for the chosen value of K. We need to choose one of the replicates to plot. To choose the 'best' replicate we will find the minimum cross-entropy criterion for that value of K
```
bestK <- ADD YOUR VALUE HERE
ce <- cross.entropy(obj.snmf, K = bestK)
best.run <- which.min(ce)

qmatrix = Q(obj.snmf, K =bestK, run=best.run)
barplot(t(as.matrix(qmatrix)), col=rainbow(bestK), xlab="Individual #", ylab="Ancestry", border=NA)
```
NOTE: The snmf runs are automatically saved into an snmf project directory - have a look in your working directory. The name of the snmf project file is the same name as the name of the input file with a .snmfProject extension ("genotypes.snmfProject").
An snmf project can be load in a different session.
project = load.snmfProject("genotypes.snmfProject")

# Extra Exercise for the maRsters

We can add metadata to the data frame we created. Load the relevant txt file that came with the data. If this file does not work - can you figure out why and can you figure out how to read this into R. Try solving it with R or BASH, not excel. Once you have solved this problem - we will them merge the metadata file with our data frame containing the PCs. To merge we need two data frames with a matching column - some or all of the values in the matching column must be identical in both data frames and the column must have the same header - in this case we are matching individual sample ID between the PC data frame and the metadata file. You can check the names of the dataframes using the names()function. You may need to change the column name in the metadata dataframe. Then use merge as follows

