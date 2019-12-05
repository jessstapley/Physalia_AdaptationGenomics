## LEA Exercise

During this exercise we will investigate population structure and perform genome-wide association tests using the R pakage LEA. The package includes statistical methods for estimating ancestry coefficients from large genotypic matrices and for evaluating the number of ancestral populations (snmf, pca). More information about the package can be found http://membres-timc.imag.fr/Olivier.Francois/LEA/index.htm. Also there a good tutorial is http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_1.html

The first part of this analysis is the same as the exercise that what we did for the 'PopulationGeneticStructure_R.md'

Open R and set up the environment.
```
library(LEA)
library(RColorBrewer)
dir.create("LEA_analyses")
setwd("/Users/jessicastapley/Dropbox/WorkshopAdaptationGenomics/CourseProgram/JS_E7_GEA/LEA_analyses")

```
You have access to .lfmm and .env files for this exercise. LEA can read multiple file types (e.g. STRUCTURE using ``` struct2geno()```. You will need to change file.lfmm to the filename I have given you for the analysis.

We will begin by performing a PCA to get an idea about the populaiton strcuture in the data. Looking at the number of principle components that explain varaition in the genomic data we can get an idea about the likely number of ancestral populations. This can help us choose a range of K for the next anaysis step.

```
# pca1R <- pca("file.lfmm", scale =TRUE)
# pca1R <- pca(paste0(input.path,".geno"), scale =TRUE)

# Plot eigenvalues.
plot(pca1R, lwd=5, col="blue", cex = .7, xlab=("Factors"), ylab="Eigenvalues", xlim=c(0,10), typ="b")

# PC1-PC2 plot.
plot(pca1R$projections)
# PC3-PC4 plot.
plot(pca1R$projections[,3:4])
```
Q1. How many ancestral populations do you think there are in your data?

Next we will investigate the population structure more throughly. The ```snmf``` function estimates ancestry coefficients similar to the commonly used programs (STRUCTURE  and ADMIXTURE). The function estimates entropy criterion that evaluates the quality of fit of the statistical model to the data using a cross-validation technique and can help you choose the number of  population clusters there are in the data. We are performing ten replicates``` rep=10```. This fuction can take some time but in comparison to running STRUCTURE it is much faster. We explore a range of values of Ks that span the true number of ancestral populations. You will need to choose a range of Ks based on your answer to Q1. E.g. if the K =3, then you could choose K=1:6.

```
obj.snmf <- snmf("file.lfmm", K=1:6, rep=5, entropy=T, ploidy =2, project ="new", iterations=10000)
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
summary(obj.snmf)
```
The plot shows the cross-entropy criterion plotted against different values of K. The lowest value of cross-entropy criterion indicates the most likely number of ancestral populations in this data. Determining the correct value of K can be difficult under some demographies. If the line is declining this would indicate finer popuation structure possibly due to isolation by distance. The summary option shows you the min, mean and max of the cross entropy values for each K. 

Q2. How many ancetral populations are evident in your data?

We will now plot the ancestry coefficients for each individual in a bar plot for the chosen value of K. We need to choose one of the replicates to plot. To choose the 'best' replicate we will find the minimum cross-entropy criterion for that value of K.
```
bestK <- 4
ce <- cross.entropy(obj.snmf, K = bestK)
best.run <- which.min(ce)

qmatrix = Q(obj.snmf, K =bestK, run=best.run)
barplot(t(as.matrix(qmatrix)), col=rainbow(bestK), xlab="Individual #", ylab="Ancestry", border=NA)
```
To optimize this plot we need to sort the individuals into their populations. We skip this for now.

NOTE: The snmf runs are automatically saved into an snmf project directory - have a look in your working directory. The name of the snmf project file is the same name as the name of the input file with a .snmfProject extension ("genotypes.snmfProject").
An snmf project can be load in a different session using 
```project = load.snmfProject("genotypes.snmfProject")```

# Genome-wide Association Scan
This LEA package can do a genome-wide association scan using Latent Factor Mixed Models (LFMM) - to find associations with the environment (ecology) or phenotype. With this approach we can correct for confounding effects of population structure due to shared demographic history or background genetic variation and find association between allele frequencies and an ecological or phenotype predictors. With large data sets this step can be computationally demanding and you might need to do your analysis on a server.

This takes some time - so for this exercise we start this running and come back to it after the lecture
```
obj.lfmm = lfmm("file.lfmm", "file.env", K = bestK, rep = 5, project="new")
```
Here we use our estimate of the best K, but is not necessarily the best approach. The value of K or latent factors to include in the model is used to help control for population structure and reduce the false discovery rate (FDR) while maintaining some power to identify outliers. So you should use several runs for each value of K and combine p-values. The manual suggests a good approach is to use the median z-scores of 5-10 runs and re-adjust the p-values afterwards to increase the power of lfmm tests (see the manual for more details). The manual also recommends using a large number of cycles (e.g., -i 6000) and the burnin period should be set to at least one-half of the total number of cycles (-b 3000). If your data set is on the small side (eg, a few hundreds of individuals, a few thousands of loci), they recommend increasing the burnin period and the total number of cycles in this situation.

When this loaded we can will adjust the z-scores and then calculate pvalues. We will use the genomic inflation factor (lambda, λ) to adjust the z-scores. Population stratification can increase the FDR in association tests. To control for this we can calculate an inflation factor, this is an estiamte of how inflated the statistic is by population stratificaiton. We estimate this for loci not linked to the trait. The method rests upon the assumptions that the inflation factor λ is constant, which means that the loci should have roughly equal mutation rates, should not be under different selection in the two populations, and the amount of Hardy-Weinberg disequilibrium measured in Wright’s coefficient of inbreeding F should not differ between the different loci. 

```
p = lfmm.pvalues(obj.lfmm, K = bestK)
pvalues = p$pvalues

zs = z.scores(obj.lfmm, K = bestK)
zs.median = apply(zs, MARGIN = 1, median)

lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

```
Q3. What is the value of lambda?

If lambda is less than or equal to 1, then no correction is needed. If you correct using a lambda less than one you can increase the FDR. We then use lambda to adjust the pvalues and plot a histogram of the pvalues. The pvalue distribution should look flat with a peak close to zero. 
```
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
hist(adj.p.values, col = "red")
```

If you are happy that you have corrected the z-scores effectively, then we can do the final setp whihc is to control for multiple testing. There are a couple of options, one uses Benjamini-Hochberg (BH) procedure.

```
## L = number of loci
L = 500
#fdr level q
q = 0.1
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
```
The other approach is to use the Bioconductor package 'qvalue' (https://github.com/StoreyLab/qvalue) and correct for multiple testing using Storey's algorithm.  

```
library(qvalue)
plot(qvalue(adj.p.values))
candidates.qv = which(qvalue(adj.p.values, fdr = .1)$signif)
```

Q4. How many candiadtes do you have using these two methods. What does this tell you about the methods?
