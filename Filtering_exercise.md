# Filtering variants in the VCF file 

Sequencing data, alignments and variant calling is a messy business. The variant call format (vcf) file that you end up with at the end of the variant calling pipeline will contain many false SNPs and INDELs. In this exercise we will consider some useful SNP filtering steps to get your varaint call set ready for further analysis. Here we use data that was produced using the GATK Variant Calling Pipeline https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145. 

My filtering recommendations are partly based on my own experience but also following great advice from others - standing on the shoulders of giants - in particular, GATK https://software.broadinstitute.org/gatk/documentation/article?id=11069, Heng Li paper https://academic.oup.com/bioinformatics/article/30/20/2843/2422145 and Jon Puritz's SNP Filtering tutorial https://www.ddocent.com/filtering/. I recommend your take a look at these webpages.

You have been assigned a vcf file (denoted as file1.vcf). When you run these commands you will need to change file.vcf to the FILENAME I gave you, or alterantively you can create a symbolic link to the file using the following command.
```
ln -s  FILENAME.vcf file1.vcf
```
Look at the vcf file ```less file1.vcf```
To quit out of 'less' press 'q'

Familiarizse yourself with the file format, see https://software.broadinstitute.org/gatk/documentation/article?id=11005 for more details. 
What FORMAT and INFO fields are included in the data. 
Is the data haploid or diploid? 

# Filtering 
We begin with some basic filtering and remove individuals that did not sequence well, i.e. have a lot of missing data. We use vcftools to create a new file with a count of the numebr loci missing for each indvidual.
```
vcftools --vcf file1.vcf --missing-indv --out file1
```
This creates a file with extension .imiss. ```less``` this file and take a look. The 5th column has the proportion of missing genotypes. We want to remove individuals with more than 99% missing. First we will use ``` awk``` to make a new file containing the names of the individuals that are missing more than 99% genotype calls.
```
awk '$5 > 0.99' file1.imiss | cut -f1 > file1_lowDP.indv
```
In the next step will combine a few filtering arguments to create a new vcf file.  
First we remove the individuals will high missing data, then we remove indels and then we remove loci that are not polymorphic.
```
vcftools --vcf file1.vcf --remove file1_lowDP.indv  --remove-indels --maf 0.0001 --recode --recode-INFO-all --out file1_f1
```
--remove file1_lowDP.indv this remove the individuals that did not sequence well 
--remove-indels pretty self explanatory
--maf 0.0001 setting the Minor Allele Frequency greater than 0.0001, which effectively remove any monomorphic sites
--recode flag tells the program to write a new vcf file
--recode-INFO-all keeps all the INFO flags from the old vcf file in the new one 
--out designates the name of the output

How many markers and individuals are left?

# Plot your data!!!
I think it always a good idea to plot and expore your data. Look at various statistics and think about possible sources of error. I would suggest choosing some variants sites (SNPs and INDELs) and look at the alignments using something like IGV. Get a feeling for what is a good variant call and what is not. Mapping quality and depth are good filters but the values you use will be different for different data sets. 

Now lets can get some summary information about the loci using bcftools and vcftools.
```
vcftools --vcf file1_f1.recode.vcf --missing-site --out file1_f1
vcftools --vcf file1_f1.recode.vcf --site-mean-depth --out file1_f1
vcftools --vcf file1_f1.recode.vcf --SNPdensity 10 --out file1_f1
vcftools --vcf file1_f1.recode.vcf --freq2 --out file1_f1
bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/MQ\t%INFO/DP\t%INFO/FS\t%INFO/QD\t%INFO/SOR\n' file1_f1.recode.vcf > file1_f1.loci
```
These vcftools commands extract various metrics that will be useful in filtering. 
--site-mean-depth
Generates a file containing the mean depth per site averaged across all individuals. This output file has the suffix ".ldepth.mean".
--freq2
Outputs the allele frequency for each site in a file with the suffix ".frq". The second option is used to suppress output of any information about the alleles.
--missing-site
Generates a file reporting the missingness on a per-site basis. The file has the suffix ".lmiss".
--SNPdensity <integer>
Calculates the number and density of SNPs in bins of size defined by this option. The resulting output file has the suffix ".snpden".

We also use bcftools to extract the Loci information - bcftools is generally much faster than vcftools. 
MQ "RMS Mapping Quality"
DP "Approximate read depth; some reads may have been filtered"
FS "Phred-scaled p-value using Fisher's exact test to detect strand bias"
QD "Variant Confidence/Quality by Depth"
SOR "Symmetric Odds Ratio of 2x2 contingency table to detect strand bias" 

Look at the files you created. HINT use ```ls -haltr``` to see the most recent files in a directory.
Do they have column headers or empty rows? This is important when reading it into R.

We need to edit the file1_f1.loci files before we can read it into R becasue some columns have only "." and this will cause R to read them as a factor and not a numnber.
```
sed 's:[\t]\.:\tNA:g' file1_f1.loci  > file1_f1.newloci
```

# Lets look at the Loci
Below are two ways to look at the loci statistics in R. The first uses the files we just made, then other uses the package ```vcfR ``` which has some nice plotting functions. First lets do it manually, then we will take a look at ```vcfR```

```
R --vanilla
lmiss1 <- read.table("file1_f1.lmiss", header=TRUE)
l1 <- read.table("file1_f1.newloci")
dp1 <- read.table("file1_f1.ldepth.mean", header=TRUE)
den1 <- read.table("file1_f1.snpden", header=TRUE)
```
Look at these files, pay attention to the column headers, we need these later. We also need to add column header to 'l1'
```
names(l1) <- c("chr", "pos", "qual", "mq", "dp", "fs", "qd", "sor")
```
We will want to look at loci missingness, lets begin there.
```
hist(lmiss1$F_MISS)
```
If you are not using the terminal in interactive mode then you need to write this file to a pdf as follows. For each plot be sure to give it a name that is descriptive. This plot will be saved to your working directory -  if your are not sure what/where your working directory is type getwd()
```
pdf("Loci_miss_Hist.pdf")
hist(lmiss1$F_MISS)
dev.off()

```
What is the maximum missingness? 

Depth filtering is a key step. Loci with very high depth (relative the to rest) are likely to be mapping errors. In this next step we plot depth across the chromosome.
```
plot(dp~pos, l1)
```
Is there any obvious outliers and are they clumped?

We can also look at the relationship betwen depth and varaince in depth. Again I think it useful to look at outliers, regions with high variance in depth may repesent mapping errors. 

```
plot(VAR_DEPTH~MEAN_DEPTH, dp1)
```
We can also limit the range of the axis to look at the spread of the data at lower values of depth and quality 
```
plot(VAR_DEPTH~MEAN_DEPTH, dp1, ylim=c(0,500), xlim=c(0,50))
```
We expect a relationship between the mean and the variance and that outliers in this plot may be problematic. 

Now lets look at the relationship betwen quality and depth. High coverage can lead to inflated locus quality scores. Heng Li showed that the relationship between depth and quality differed between real and spurious variant calls. Considering variants with a mean depth greater than the mean depth plus 2-3 times the square root of mean depth, the for quality score should be twice as large for the real variants. Jon Puritz considered this too conservative for reduced representation libraries. 

```
plot(qual~dp, l1)
plot(mq~dp, l1)
plot(qd~dp, l1)
```
What is the mean and median depth in 'l1'? How does that compare to mean depth in the 'dp1'? 

We can also limit the range of the axis to look at the spread of the data at lower values of depth and quality. 
```
plot(qual~dp, l1, ylim=c(0,1000), xlim=c(0,1000))
```

We will want to remove variants in high depth as these are likely due to copy number variation or other larger structural events.

Now take a look at the SNP density file. What is the distribution of variants per KB, or number of SNPs ini 10 bp window. Are any of these values unusually high?

```
hist(den1$SNP_COUNT)
hist(den1$VARIANTS.KB)
plot(SNP_COUNT~BIN_START, den1)

```
I dont think there is too much to worry about, but depending on your data you may want to consider this as a possible filtering step.

Lets combine the data frames to make it a little easier to campare some metrics.
If the the order of markers is the same we can simply 'cbind' these two dataframes.

```
ldp1 <- cbind(l1, dp1)
```
We can test if the positions are the same using an ifelse statement, if the POS and pos are not the same this will produce a 1. If we sum this it will give us a total number of differences
```
sum(ifelse(ldp1$pos == ldp1$POS, 0,1))
```
This is 0 so we can continue.

```
plot(MEAN_DEPTH~dp, ldp1)
plot(qd~MEAN_DEPTH, ldp1)
range(ldp1$qual, na.rm=TRUE)

```
What is the minimum quality value for this data? Do we need to filter for minimum quality?

# Enough exploring, lets make some hard filtering decisions.

There is not one ring to rule them all - and by ring I mean hard filters. GATK guidelines are a good start. After implimenting these filters I would agian plot the data and look at aligments in IGV. 
E.g. QD < 2.0 QUAL < 30.0 SOR > 3.0 FS > 60.0 MQ < 40.0 MQRankSum < -12.5 ReadPosRankSum < -8.0
If you use GATK pipeline then use their functions (https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php)

I would begin with a depth filter - remove very low and very high depth. The choice of depth depends on your coverage. What is the depth range? It appears that the minimum depth is quite high in this dataset so we dont need to use a minimum depth filter. If you did I would suggest minimum depth of 3.

For our max depth filter lets look at the quantiles and remove outliers.
```
quantile(ldp1$dp, probs=seq(0,1,0.1))
```
Lets filter the loci with depth greater than the 90% quantile
```
ldp1.fdp <- subset(ldp1, ldp1$dp<quantile(ldp1$dp, probs=0.9))
dim(ldp1)
dim(ldp1.fdp)
```
How many loci did we loose?
Would we get the same results if we filter with mean depth?

```
ldp1.fdp <- subset(ldp1, ldp1$MEAN_DEPTH<quantile(ldp1$MEAN_DEPTH, probs=0.9))
dim(ldp1)
dim(ldp1.fdp)
```
Is there much difference? 

Now lets replot the data.
```
plot(qual~dp, ldp1.fdp)
plot(qd~dp, ldp1.fdp)
plot(mq~dp, ldp1.fdp)
plot(VAR_DEPTH~MEAN_DEPTH, ldp1.fdp)

```
Do you see any outliers in these plots? 

Now we remove loci based on the GATK suggestions.
```
ldp1.filt <- subset(ldp1.fdp, ldp1.fdp$qd>2 & ldp1.fdp$sor<3 & ldp1.fdp$fs<60 & ldp1.fdp$mq>40)
dim(ldp1.filt)
```
How many loci are left?

Now look at the plots again.

```
par(mfrow=c(1,2))
plot(qual~dp, ldp1.filt, main="Depth + GATK Filtered")
plot(qual~dp, ldp1.fdp, main="Depth Filtered")

par(mfrow=c(1,2))
plot(qd~dp, ldp1.filt, main="Depth + GATK Filtered")
plot(qd~dp, ldp1.fdp, main="Depth Filtered")

par(mfrow=c(1,2))
plot(mq~dp, ldp1.filt, main="Depth + GATK Filtered")
plot(mq~dp, ldp1.fdp, main="Depth Filtered")

par(mfrow=c(1,2))
plot(VAR_DEPTH~MEAN_DEPTH, ldp1.filt, main="Depth + GATK Filtered")
plot(VAR_DEPTH~MEAN_DEPTH, ldp1.fdp, main="Depth Filtered")

```
If there are any loci that have high variance for their mean depth, I would remove these using something like this.

```
ldp1.filt$vmdp <- ldp1.filt$VAR_DEPTH/ldp1.filt$MEAN_DEPTH
ldp1.filt1 <- subset(ldp1.filt, ldp1.filt$vmdp<quantile(ldp1.filt$vmdp, prob=0.99))

dim(ldp1.filt1)

par(mfrow=c(1,2))
plot(VAR_DEPTH~MEAN_DEPTH, ldp1.filt, main="Filtered")
plot(VAR_DEPTH~MEAN_DEPTH, ldp1.filt1, main="Depth Variance Filtered")

```
The last step is to write out two lists of the loci from the main filtering steps. Normally you are only interested in a single list of your best SNPS, but we create two lists here to make some comparisons with later.
```
chr.pos.dp <- ldp1.fdp[,c(1,2)]
chr.pos.f1 <- ldp1.filt1[,c(1,2)]

write.table(ldp1.fdp[,c(1,2)], file="DepthFilteredLoci.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(ldp1.filt1[,c(1,2)], file="DepthGATKFilteredLoci.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
quit()
```
We quit out of R. Now in your working directory you will find the new files you created. We will pass these to vcftools to create two filtered vcf files. We add another filtering step here to remove loci with max missing 0.1. This last additional filtering step is really just to remove loci with high missingness. Later we consider other options for removing loci with different missingness. 

```
vcftools --vcf file1_f1.recode.vcf --positions DepthFilteredLoci.txt  --max-missing 0.1 --recode --recode-INFO-all --out file1_fltdp
vcftools --vcf file1_f1.recode.vcf --positions DepthGATKFilteredLoci.txt --max-missing 0.1 --recode --recode-INFO-all --out file1_fltqatkdp

vcftools --vcf file2_f1.recode.vcf --positions DepthFilteredLoci.txt  --max-missing 0.1 --recode --recode-INFO-all --out file2_fltdp
vcftools --vcf file2_f1.recode.vcf --positions DepthGATKFilteredLoci.txt --max-missing 0.1 --recode --recode-INFO-all --out file2_fltqatkdp
```
How many loci were retained?

After hard filtering it is a good idea to create summary files to double check everythig has worked as planned and you can use these files later when describing the data.
```
vcftools --vcf file1_fltdp.recode.vcf --missing-indv --out file1_fltdp
vcftools --vcf file1_fltdp.recode.vcf --missing-site --out file1_fltdp
vcftools --vcf file1_fltdp.recode.vcf --site-mean-depth --out file1_fltdp
vcftools --vcf file1_fltdp.recode.vcf --SNPdensity 10 --out file1_fltdp
vcftools --vcf file1_fltdp.recode.vcf --freq2 --out file1_fltdp
bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/MQ\t%INFO/DP\t%INFO/FS\t%INFO/QD\t%INFO/SOR\n' file1_fltdp.recode.vcf > file1_fltdp.loci

vcftools --vcf file1_fltqatkdp.recode.vcf --missing-indv --out file1_fltqatkdp
vcftools --vcf file1_fltqatkdp.recode.vcf --missing-site --out file1_fltqatkdp
vcftools --vcf file1_fltqatkdp.recode.vcf --site-mean-depth --out file1_fltqatkdp
vcftools --vcf file1_fltqatkdp.recode.vcf --SNPdensity 10 --out file1_fltqatkdp
vcftools --vcf file1_fltqatkdp.recode.vcf --freq2 --out file1_fltqatkdp
bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/MQ\t%INFO/DP\t%INFO/FS\t%INFO/QD\t%INFO/SOR\n' file1_fltqatkdp.recode.vcf > file1_fltqatkdp.loci
```

Lets plot again to check everything worked.

```
R --vanilla
imiss <- read.table("file1_fltqatkdp.imiss", header=TRUE)
lmiss <- read.table("file1_fltqatkdp.lmiss", header=TRUE)
dp <- read.table("file1_fltqatkdp.ldepth.mean", header=TRUE)
den <- read.table("file1_fltqatkdp.snpden", header=TRUE)
l1 <- read.table("file1_fltdp.loci")

hist(lmiss$F_MISS)
hist(imiss$F_MISS)
```

What is the maximum per individual missing rate?
What is the maximum per site missing rate? 

You can do some other exploratory plots if you have time.

```
quit()
```

Now we have a filtered set of loci that can be used for further analysis. Note for downstream analyses you will be required to do additional filtering, for example different levels of missingness per individual and per loci, as well exlcuding markers with low minor allele frequencies.

# Population specific Loci filtering 
If you know a priori that you have genetically distinct populations and your analysis will rely on comparisons across populations, then it can be useful to look at the loci missingness across populations and implement population specific filtering steps. For example, if a loci has high missingness in a single population, then it can be better to exlcude this loci when comparing between populations. If only one or two individuals in a population have reads covering this site, then this can skew the allele frequency in a population.

We do not have population information for this data set so we will not be doing this filtering step - I just put this infromation here so you have it.

Here is an example of how to use a population specific loci missingness filter. The proportion of missingness (cutoff) that you choose at this stage is dependent on the number of individuals in the population.

We can use vcftools to subset the vcf file to a single population and then calculate the site missingness. The id.pop file is a list of indivduals in a sigle population.  
```
vcftools --vcf file.vcf --keep id.pop1 --missing-site --out file_pop1
vcftools --vcf file.vcf --keep id.pop2 --missing-site --out file_pop2
vcftools --vcf file.vcf --keep id.pop3 --missing-site --out file_pop3
```
This creates three .imiss files, we then combine these and select a cutoff to create a list of loci that do not perform well cross your populations - here we choose 0.2. 
```
cat *.lmiss | awk '!/CHR/' | awk '$6 > 0.2' | cut -f1,2 >> badloci 
vcftools --vcf file.vcf --exclude-positions  badloci --recode --recode-INFO-all --out file_popG80
```

# Visualizing and Filtering with vcfR
vcfR package provides some nice features for working with vcf files in R and I recommend you have a look at the documentation created by Brian J. Knaus and Niklaus J. Grünwald https://knausb.github.io/vcfR_documentation/index.html
The infromation below is based on their Filtering tutorial https://knausb.github.io/vcfR_documentation/filtering_data.html

A nice feature of the package is that it can plot the SNP information alongside annotation information. It might seem like this package makes what we did before redundant, but I think it is still important to do the exploratory analysis 'manually' to better understand variant quality. Also vcfR can only handle chromosome level data - so if you dont have a good reference genome you may encounter some problems.

If you have a vcf file and gff files for the entire genome you will need to subset the data into chromosomes (see https://knausb.github.io/vcfR_documentation/subset_data_to_1chrom.html)

First we load the package ```vcf``` and the data.
```
vcf_file <-  "file.vcf"
vcf <- read.vcfR( vcf_file, verbose = FALSE)
```
Then we create a chromosome object and plot basic varaint statistics across the chromosome.
```
chrom <- create.chromR(name='groupI', vcf=vcf)
plot(chrom)

```
Check out that plot!

You can use the masker funtion to filter SNPs. Below are recommendations from the tutorial, how do this compare to your data?
```
chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = FALSE)
plot(chrom)
```
If your happy with these filtering decision you can write out a filtered vcf file
```
write.vcf(chrom, file="good_variants.vcf.gz", mask=TRUE)
```
