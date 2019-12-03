# Investigating Population Genetic Struture in R
In this exercise we invesitgte population genetic strcuture using two approaches - a distance based method (PCA) and a model based method (entropy criterion).

First we is create a new directory for the R project - call it 'PopGen'. In that directory create one subdirectory called data. Copy the data file to that data directory.

Data files
Brown: SimInversion

Green: Chicken

Red: Cichlid

Blue: Salmon

Open R studio, create a new project using this existing directory that you just created (PopGen). Once R project is et set you have two options. if you are confortable in R make an R markdown file. If you are a beginner in R then just make a new script.

For R script - you can copy and paste the chunks of code from this github into the R script. In the Rscript you can add comments and notes - use the '#' before a line of comments. This will help you see what is code and what are comments.

For R Markdown - in R open a new R Markdown file. Have a look at the default file format. You can add chuncks of code to the sections between "```"  "```". You can replace the code in all sections EXCEPT this bit - DO NOT OVERWRITE this bit.
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
Get familiar with this sytax and organisation of the markdown documents. Save this document (PopGenStructure) to the directory PopGen and now you can edit this document by copying the information from this github tutporial - copy the code from here into the markdown document in chunks. To run a chunk of code in the markdwon document press the small green arrow.

# PCA analysis from vcf file
Now the R project is set up lets begin.
Read in your libraries.

```
library(LEA)
library(RColorBrewer)

```
You will need to change the name of your input.path and out.path in the code below to the file that is you data directory (HINT look at the 'files'). First we convert the vcf to geno format for this analysis
```
input.path <-  "data/FILENAME
vcf2geno(paste0(input.path,".vcf"),paste0(input.path,".geno"))
```
Then we performa PCA and plot the data
```
pca1R <- pca(paste0(input.path,".geno"), scale =TRUE)
par(mfrow=c(2,2))

# Plot eigenvalues.
plot(pca1R, lwd=5, col="blue", cex = .7, xlab=("Factors"), ylab="Eigenvalues")

# PC1-PC2 plot.
plot(pca1R$projections)
# PC3-PC4 plot.
plot(pca1R$projections[,3:4])

# Plot standard deviations.
plot(pca1R$sdev)
```
Looking at the number of principle components that explain variation in the genomic data we can get an idea about the likely number of ancestral populations. This  can help us choose a range of K for the next anaysis step.

The number of ”significant” components can be evaluated using graphical methods based on the screeplot. The knee in the screeplot indicates that number of genetic clusters). 

```
screeplot(pca1R)
plot(pca1R)# plots the eigenvalues

```
Next we use calculate Tracy-Widom test for each eigenvalue to investigate the number of genetic clusters in the data. We display the p-values for the Tracy-Widom test for 1-10 eigenvalues. The 'knee' in the plot indicates the number of significant components (K) in the data. The number of genetic clusters is K+1.

```
tw <- tracy.widom(pca1R)
tw$pvalues[1:10]
plot(tw$percentage, pch=19, typ="b") 
```
Q1. How many populations do you think there are in your data?

# Looking at population structure 

Next we will investigate the population structure more throughly. The ```snmf``` function estimates ancestry coefficients similar to the commonly used programs (STRUCTURE  and ADMIXTURE). First we convert the data again.

```
geno2lfmm(paste0(out.path,".geno"))

```
Look at the files in the data directory

Now do the population clustering. The function estimates entropy criterion that evaluates the quality of fit of the statistical model to the data using a cross-validation technique and can help you choose the number of population clusters there are in the data. We are performing ten replicates``` rep=10```. This fuction can take some time but in comparison to running STRUCTURE it is much faster. We explore a range of values of Ks that span the true number of ancestral populations. You will need to choose a range of Ks based on your answer to Q1. E.g. if the K =3, then you could choose K=1:6.

This function takes to long to run so we will load premade data. You need to look in your files and find out waht the .snmf project is called and replace FILENAME with the name of your file.

```
# this is how you woud run the code
# obj.snmf <- snmf(paste0(input.path,".lfmm"), K=1:6, rep=10, entropy=T, ploidy =2, project ="new", iterations=10000)
# this is how re reload a prepared object
obj.snmf = load.lfmmProject("FILENAME.lfmmProject")

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

How could we make these plots better, what information are we missing?

