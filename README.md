# Functional enrichment for QTL variants 

This script performs functional enrichment of QTL variants.

## Usage example 

The functional enrichment is performed in two steps. The first step creates a file containing for each variant inside a vcf, the MAF and the closest upstream and downstream feature (e.g gene). Second, it will read through all QTL variants and for each of them randomly pick variants (from the file created previously) that match MAF (±2%) and distance to TSS (±2.500) and that are not eQTLs. This allows for the generation of a null distribution. Then, it uses a Fisher exact test to perform enrichment over the null of the QTL variants.

## Step by step example

### Step1: Create null file 

Before creating the null file you also need to provide a list of variants that are nominally significant eQTLs. These can be provided by using QTLtools cis --nominal 0.05. Then you need to only keep the var_id. It should be the 8th column.


```bash
fenrich null \
  --vcf genotypes.chr22.vcf.gz \
  --bed genes.50percent.chr22.bed.gz \
  --qtl nominal_only_significant_snps.txt.gz \
  --out null_variantFile.txt.gz
```



### Step2: Run enrichment analysis 

**IMPORTANT: the eQTL results should be in QTLtools cis format. Please check QTLtools for the exact format to use.**


```bash
fenrich enrich \
  --null null_variantFile.txt.gz \
  --phen TFs.encode.bed.gz \
  --qtl results.permute.genes.significant.txt.gz \
  --out qtlTFenrichment.txt \
  --random_var 100 \
  --window_size 2500 \
  --maf_window 0.02

```

The output of fenrich enrich mode looks like this: 

```
#A       B       C       D
#59      22      214     531
```
A -> Number of eQTLs overlapping with the peaks. 

B -> Number of null variants overlapping with peaks. 

C -> Number of eQTls not overlapping with peaks. 

D -> Number of null variants not overlapping with peaks. 

Next, you can load the results file into R, create a matrix and run a fisher exact test as shown below. 

```{r}
D = read.table("qtlTFenrichment.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
mat = matrix(D,ncol=2, byrow=TRUE)
fisher.test(mat)
```

