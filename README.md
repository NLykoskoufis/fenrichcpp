# Functional enrichment for QTL variants 

This script performs functional enrichment of QTL variants.

## Usage example 

The functional enrichment is performed in two steps. The first step creates a file containing for each variant inside a vcf, the MAF and the closest upstream and downstream feature (e.g gene). Second, it will read through all QTL variants and for each of them randomly pick variants (from the file created previously) that match MAF (±2%) and distance to TSS (±2.500) and that are not eQTLs. This allows for the generation of a null distribution. Then, it uses a Fisher exact test to perform enrichment over the null of the QTL variants.

## Step by step example

### Step1: Overlap using Bedtools

First, we need to overlap our variants with the ChIP-seq data we want to perform enrichment.

```{bash}
bedtools intersect -a variants.bed.gz -b ChIP_seq_data.bed.gz -wa -wb > overlaped_elements.txt
```
The output file should look like this

```
1       762272  762273  rs3115849       G/A     +       1       762183  762956  CTCF
1       762588  762589  rs71507461      G/C     +       1       762183  762956  CTCF
1       762591  762592  rs71507462      C/G     +       1       762183  762956  CTCF
1       762600  762601  rs71507463      T/C     +       1       762183  762956  CTCF
1       762631  762632  rs61768173      T/A     +       1       762183  762956  CTCF
1       805555  805556  rs72631880      T/A     +       1       804870  805987  CTCF
1       839872  839873  rs192553893     C/T     +       1       839340  840863  CTCF
```

### Step2: Create null file 

Before creating the null file you also need to provide a list of variants that are nominally significant eQTLs. These can be provide by using QTLtoos cis --nominal 1. Then you need to only keep the significant snps, like so

```{bash}
zcat nominal1_chrALL.txt.gz | awk '{if($12 <0.05) {print $8}}' > nominal_only_significant_snps.txt
```

```{bash}
fenrich null \
    --vcf <vcf file> \
    --bed <bed file> \
    --qtl <QTLs> \
    --out <output file>
```

In order to speed the process you can also perform this by chromosome

```{bash}
fenrich null \
    --vcf <vcf file> \
    --bed <bed file> \
    --qtl <QTLs> \
    --out <output file> \
    --region <chromosome>
```

### Step3: Run enrichment analysis 

By default the window size is set to ± 2'500, the maf window to ±2% and the number of random variants to 10. You can always change these parameters if needed.

```{bash}
fenrich enrich \
    --nul <null file > \
    --qtl <qtls to be enriched> \
    --phen <variant overlaped with peaks file> \
    --out <out file> 
```


