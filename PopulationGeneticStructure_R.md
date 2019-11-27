# Investigating Population Genetic Struture in R
This exercise we invesitgte population genetic strcuture using two approaches - a distabce based method (PCA) and a model based method (entropy criterion).

This package can read in multiple file formats, including vcf, plink ped, STRUCTURE. Here I provide an example code for vcf and lfmm file formats.

```
library(LEA)
library(RColorBrewer)

dir.create("LEA_analyses")
setwd("LEA_analyses")
```
You will need to change the name of your file in the code below.
```
# with .lfmm file
pca1R <- pca("file.lfmm", scale =TRUE)

# with .vcf file you convert to geno first
vcf_file <- "file.vcf"
vcf2geno(vcf_file, output.file = "file.geno", force = TRUE)      
```
We will begin by performing a PCA to get an idea about the populaiton strcuture in the data. Looking at the number of principle components that explain variation in the genomic data we can get an idea about the likely number of ancestral popolations. This 
can help us choose a range of K for the next anaysis step.

```
pca1R <- pca("file.lfmm", scale =TRUE)
# or for vcf 
pca1R <- pca("file.geno", scale =TRUE)
tw <- tracy.widom(pca1R)
tw$pvalues[1:10]
plot(tw$percentage, pch=19, typ="b") 

```
Here we perform a pca and then we use the Tracy-Widom test for each eigenvalue to investigate the number of genetic clusters in the data. We display the p-values for the Tracy-Widom test for 1-10 eigenvalues. The 'knee' in the plot indicates the number of significant components (K) in the data. The number of genetic clusters is K+1.

Q1. How many ancestral populations do you think there are in your data?

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

qmatrix = Q(obj.snmf, K =4, run=best)
barplot(t(as.matrix(qmatrix)), col=rainbow(bestK), xlab="Individual #", ylab="Ancestry", border=NA)
```
To optimize this plot we need to sort the individuals into their populations. We skip this for now.

NOTE: The snmf runs are automatically saved into an snmf project directory - have a look in your working directory. The name of the snmf project file is the same name as the name of the input file with a .snmfProject extension ("genotypes.snmfProject").
An snmf project can be load in a different session.
project = load.snmfProject("genotypes.snmfProject")
