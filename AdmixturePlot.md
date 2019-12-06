
# Creating Plot in R of Admixture output 
Admixture was run using .bed files. These are created using plink. The admixture command line was

```
for K in 1 2 3 4 5 6 7 8 9 10; \
do ../../Admixture/admixture --cv=10 cichlid.bed $K | tee log${K}.out; 
done
```

We run multiple replicates of admixture across multiple values of K and use cross validation (CV) to  detrmine the most likely number of genetic clusters (K). 
The K with the lowest cross-validation error compared to others is the most likely number of ancestral clusters

We can look at the CV values in the log file by extracting the relevant information from the log files using grep. In the directory that contains the output from admixture type
```
grep -h CV log*.out
```

E.g. 
grep -h CV log*.out

CV error (K=5): 0.59575

CV error (K=6): 0.59378

CV error (K=7): 0.59293

CV error (K=8): 0.59239 ***

CV error (K=9): 0.59353

CV error (K=10): 0.59365

The lowest value is the most likely K, in this e.g. K=8

# Plotting the admixture results

The admixture results are in the same format as the STRUCTURE and snmf outputs. The .Q file is an ancestry matrix, each indvidual is a row and there are K number of columns. The values are the proportion of inferred ancestry for that genetic cluster.

We will use R to plot this data. Open R and load at least 3 x .Q files, for the best value of K (lowest CV) and the file for K-1 and K+1. 
The following was written assuming K=8 is the K with lowest CV, the output files are in a folder in your projet directory called 'data' and your files are called e.g. file1.7.Q. You will need to change the filenames but also think about how you need to change the number of colours on the plot (e.g. rainbow(8))

```
tbl1 = read.table("data/file1.7.Q")
tbl2 = read.table("data/file2.8.Q")
tbl3 = read.table("data/file3.9.Q")


par(mfrow=c(3,1))
barplot(t(as.matrix(tbl1)), col=rainbow(7), xlab="Individual #", ylab="Ancestry", border=NA)
barplot(t(as.matrix(tbl2)), col=rainbow(8), xlab="Individual #", ylab="Ancestry", border=NA)
barplot(t(as.matrix(tbl3)), col=rainbow(9), xlab="Individual #", ylab="Ancestry", border=NA)

```

