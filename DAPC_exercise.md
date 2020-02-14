# Discriminant analysis of principal components (DAPC)

We will perform DAPC using the R pacakge ```Adegenet```. The vignette can also be opened within R by executing ```adegenetTutorial("dapc")```.
Copy the relevant file to your R project directory. This appraoch is interactive - you will be promted to enter values in the R console during the analysis. The choise of the values is dependent on the plots that the functions produce.

Open R and set up a R PROJECT for this exercise and load the relevant libraries.

```
library(vcfR)
library(adegenet)

```
Read in the vcf file, remember you will need to replace "file.vcf" with the name of your vcf file. 
```
vcf_file <- "file.vcf"
gld <- vcfR2genind(vcf_file)
```
The first step is to find clusters in the data, we set a maximum number of clusters to be some value exceedng what we expect is present in the data.
```
grp.ld <-  find.clusters(gld, max.n.clust=30)
```
In this step the function performs a PCA and creates a plot of the cumulative variance explained by increaseing values of PC retained. At this point you are asked to choose the number of PC to retain. Keep all the PCs at this stage, so enter the maximum value on the x-axis. 

Then the function plots a graph of the Bayesian Information Criteron (BIC) for different values of K. The elbow in the curve and the lowest BIC is the most likely estimate of the number of clusters in the dataset. You are again promted to enter a value, this time a value of the likely number of clusters is needed.

In the next step we plot a DAPC. Again you are aksed to provide estimats of the number of PCs to retain - this time choose a lower value - where the curve plateaus (i.e. where more PC do not explain more variation) and choose the value of number of clusters 
```
dapc.ld <- dapc(gld, grp.ld$grp, n.da = NULL)
scatter(dapc.ld)

```

