# Functional enrichment for QTL variants 

This script performs functional enrichment of QTL variants.

## Usage example 

The functional enrichment is performed in two steps. The first step creates a file containing for each variant inside a vcf, the MAF and the closest upstream and downstream feature (e.g gene). Second, it will read through all QTL variants and for each of them randomly pick variants (from the file created previously) that match MAF (±2%) and distance to TSS (±2.500) and that are not eQTLs. This allows for the generation of a null distribution. Then, it uses a Fisher exact test to perform enrichment over the null of the QTL variants.

## Step by step example

### Step1: Check that phenotype file is sorted.

First, we need to check that the phenotype file is sorted. If not sorted then please sort using either unix-sort or bedtools sort. The file should be in bed format without the header.

```bash
bedtools sort -i H3K4me3_data.bed > H3K4me3_data_sorted.bed
# or 
sort -V -k1 -k2,2n -k3,3n H3K4me3_data.bed > H3K4me3_data_sorted.bed
```


### Step2: Create null file 

Before creating the null file you also need to provide a list of variants that are nominally significant eQTLs. These can be provided by using QTLtools cis --nominal 0.05. Then you need to only keep the var_id. It should be the 8th column
```bash
zcat nominal005_chrALL.txt.gz | cut -d" " -f8 | sort | uniq > nominal_only_significant_snps.txt
```

```bash
fenrich null \
    --vcf <vcf file> \
    --bed <bed file> \
    --qtl <QTLs> \
    --out <output file>
```

In order to speed the process you can also perform this by chromosome

```bash
fenrich null \
    --vcf <vcf file> \
    --bed <bed file> \
    --qtl <QTLs> \
    --out <output file> \
    --region <chromosome>
```

--vcf | vcf file (should be the same as the one you used for eQTL discovery) !!! VCF file should contain GT information

--bed | bed file containing gene positions (should be the same bed file as the one used for eQTL discovery)

--qtl | txt file containing SNPid of nominally significant SNPs (one per line).




### Step3: Run enrichment analysis 

By default the window size is set to ± 2'500, the maf window to ±2% and the number of random variants to 10. You can always change these parameters if needed.

**IMPORTANT: the eQTL results should be in QTLtools cis format. Please check QTLtools for the exact format to use.**


```bash
fenrich enrich \
    --nul <null file > \
    --qtl <qtls to be enriched> \
    --phen <phenotypes> \
    --out <out file> 
```
--phen | the phenotypes you want to enrich your eQTL SNPs for in bed format.
