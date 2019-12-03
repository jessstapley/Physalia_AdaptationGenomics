# Multiple Methods to LD prune a vcf file 

Many analyses need the SNPs to be in approximate linkage equilibirium, thus a first step to SNP analysis is to prune the SNPs to remove loci in high linkage disequilibrium (LD). We can do this using ```plink ``` or ```bcftools```.
In both cases we begin with a vcf file. 

# Pruning with plink 
You can dowload plink from https://www.cog-genomics.org/plink/1.9/.
Then use the following command to convert the vcf file into plink file (.ped and .map) or if you add ```-make-bed```to the command it will create the binary plink files (.bed, .dim, .fam). Binary formats are preferred for larger datasets
Things to note
The ped/fam files will include a column for family ID. To set this to zero use ```--const-fid```. 
If your SNPs are positioned on chromosomes and the total number of chromoosmes is greater than 26 then you need to specify the chromosomes using ```chr-set``` for e.g. 
```
plink --vcf filename.vcf  --chr-set 38 no-xy --const-fid --keep-allele-order  --recode --out filename
```
If your SNPs are positioned on scaffolds or contigs and there are more than 95 then you can must use ``` --allow-extra-chr ```. Note that the contig/scaffold names must be alphanumerical (not just numbers and not called chr1 or chr2).
```
plink --vcf filename.vcf  --allow-extra-chr --const-fid --keep-allele-order  --recode --out filename
```
If your SNPS do not have names (only positions) you will need to give them a name, by editing the .map or .bim files. This is important for the pruning step. To prune we use ```--indep-pairwise``` option. This is followed by three numbers: a window size in variant count, a variant count to shift the window at the end of each step and the r^2 threshold. Note that here we are using the binary files and use the argument ```--bfile```
```
plink --bfile filename --chr-set 38 no-xy --indep-pairwise 50 5 0.3 --recode --out filename
```
This creates two files .prune.in and .prune.out, which are a list of SNP names that will be included/excluded from a new file. Here use the ```--extract``` option to select the markers to include.
```
plink --file filename --chr-set 29 no-xy  --keep-allele-order --extract filename.prune.in --make-bed --recode --out filename_pruned
```
Additional filtering can then be done on the LD pruned file. Multiple filtering statements can be included within a plink command, for example you could filter based on Hardy–Weinberg equilibrium (--hwe 0.001), use call rate at least 95% (--mind 0.05), filter out missing SNP above 2% (--geno 0.02), and remove low minimum allele frequency (--maf 0.05). 

```
plink --file filename_pruned --chr-set 29 no-xy --mind 0.05 --hwe 0.001 --keep-allele-order - --make-bed --recode --out filenames_filtered
```

# Pruning with bcftools
When converting vcf to plink you loose a lot of FORMAT and INFO information from the vcf file (i.e. depth, genotype quality) also the order of the alleles and ancestral allele information can be lost/altered. If these things are required in downstream analysis  - dont panic - we can use bcftools to do LD pruning. Bcftools is part of samtools, which can be download from http://www.htslib.org/download/. The LD pruning uses the prune plugin https://samtools.github.io/bcftools/howtos/plugins.html. One line and yourare done.
```
bcftools +prune -l 0.2  -w 5000 file.vcf -Ov -o file_pruned.vcf
```

# Pruning with R
Seveal packages in R can LD prune a dataset, below is an appraoch using the ```SNPRelate``` package. This requires the data in the gds format, so the first step is convert vcf to gds. 

```
library(SNPRelate)
snpgdsVCF2GDS("file.vcf", "file.gds")
```
This will output some basic information to the console, i.e. number of markers and samples. Now we can prune using the default settings. 
```
gds <- snpgdsOpen("file.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(gds)

```
This outputs information to the console, including the number of markers excluded from each chromosome. The output is a list of the SNPs to exclude from each chromosome. You can use the follwoing to look at specific chromosomes or SNPs. 
```
names(snpset)
head(snpset$chr1)
```
You can create a snp.id list as folllows
```
snp.id <- unlist(snpset)
```
