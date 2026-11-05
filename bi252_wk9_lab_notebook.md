Notebook for BI252 Lab Wk9: Genetic diversity lab
================
Chloe and Asher
2026-11-05

## Learning goals

1.  Navigate and understand the variant call format (VCF), which stores
    genotype data

2.  Calculate heterozygosity and effective population size (Ne) using
    genotype data

3.  Analyze genetic relatedness using principal components analysis
    (PCA)

## Introduction

Human populations that have migrated away from the origin of Homo
sapiens in the Addis Ababa region of Ethiopia exhibit less genetic
diversity and smaller effective population size due to founder effects.

In this lab, we will analyze genotyping data in chromosome 8 of seven
individuals from the Genome in a Bottle Project (GIAB). Below is a table
of each patient and the migration distance of their ancestry group from
Addis Ababa.

The patients below are from each ancestry group, whose expected
migratory distances from the Addis Ababa region of Ethiopia are below:

Patient Ancestry Ancestral region Migration distance

HG001 Western European Western Europe (Paris) ~5,600 km HG002 Ashkenazi
Jewish Levant (Jerusalem) ~2,600 km HG003 Ashkenazi Jewish Levant
(Jerusalem) ~2,600 km HG004 Ashkenazi Jewish Levant (Jerusalem) ~2,600
km HG005 Han Chinese Northern China (Beijing) ~8,300 km HG006 Han
Chinese Northern China (Beijing) ~8,300 km HG007 Han Chinese Northern
China (Beijing) ~8,300 km

This lab will evaluate the following prediction:

Individuals from ancestries possessing fewer migrations out of Addis
Ababa Ethiopia will exhibit greater heterozygosity and effective
population size.

We will evaluate this prediction by:

1.  Analyzing relatedness using PCA
2.  Calculating heterozygosity from the genotype data
3.  Evaluating the correlation between heterozygosity and migration
    distance

## Data

- VCFs on seattle for seven human genomes: hg001 through hg007

## Software

- seattle
- FileZilla
- IGV
- R markdown

## Instructions

### Step 1 - Create and visualize a VCF

1.  Log into seattle using ssh <username@seattle.skidmore.edu>
2.  Run the bash code below

``` bash

# create a wk 9 lab directory and replace your username
mkdir /biology/data/cnieliwocki/wk9_lab

# move to wk 9 lab directory and replace your username
cd /biology/data/cnieliwocki/wk9_lab

# load bcftools
module load anaconda3

conda activate /cm/shared/apps/anaconda3/biology/envs/bcftools

# merge your VCF with those of hg001, and hg003-hg007
# copy and run everything from merge to ".gz" below and run 
bcftools merge --force-samples \
    /biology/data/class-bi252/wk9_lab/hg001_chr8.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg002_chr8_trio.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg003_chr8.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg004_chr8.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg005_chr8.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg006_chr8.vcf.gz \
    /biology/data/class-bi252/wk9_lab/hg007_chr8.vcf.gz \
    -Oz -o hg0001_7_merged.vcf.gz
    
# index the vcf so it can be loaded into IGV
bcftools index -t hg0001_7_merged.vcf.gz
```

3.  Download data using FileZilla

- Open FileZilla
- Enter ‘seattle.skidmore.edu’ in Host
- Enter your Skidmore username in Username
- Enter your seattle password
- Type ‘22’ in Port
- Click Quickconnect
- Navigate to /biology/data/username/wk9_lab in the right hand tab
- Drag hg0001_7_merged.vcf.gz and hg0001_7_merged.vcf.gz.tbi to your
  Downloads in the left hand tab
- Move this .Rmd file and hg0001_7_merged.vcf.gz +
  hg0001_7_merged.vcf.gz.tbi into the same folder together

4.  Upload to IGV

- Open IGV
- File \> Load from file \> navigate to and select
  hg0001_7_merged.vcf.gz
- In IGV, search for chr8:50,000,035-50,000,074

**Q1 - Click on the dark blue SNV in the hg001_7_merged.vcf.gz row. What
are the reference and alternate alleles?**

``` text

Reference: C

Alternate: A
```

In the rows for each patient’s genotypes, non-reference homozygous
genotypes are colored light blue. Non-reference heterozygous genotypes
are colored dark blue. Click on each genotype.

**Q2 - What are the genotypes of patients HG001-HG007 at chromosome 8
position 50,000,055? What are their genotype quality scores?**

``` text

HG001 (GT/GQ): 790

HG002 (GT/GQ): 446

HG003 (GT/GQ): 495

HG004 (GT/GQ): 495

HG005 (GT/GQ): 231

HG006 (GT/GQ): 344

HG007 (GT/GQ): 333
```

There are two trios represented in this VCF:

HG002-4, where 3 is the father of 2 and 4 is the mother of 2

HG005-7, where 6 is the father of 5 and 7 is the mother of 5

**Q3 - At this site, is it possible that each allele in HG002 and HG005
were inherited from their parents? Explain.**

``` text

Answer: It is possible. For the former pair (hg002 and its parents) have an allele combination that allows the possibility of inheritance: the child has A/A whereas both parents have a combination of A and C. For the latter pair, it is also possible. Each patient has A/A. 
```

Navigate to chr8:20,098,227-20,098,786, where you’ll see three SNVs.
Click on one of the genotypes of HG002.

**Q4 - Does the genotype contain a `/` or `|` character between each
allele?**

``` text

Answer: Each of the SNVs have a '|'.
```

Click on one of the genotypes of HG005.

**Q5 - Does the genotype contain a `/` or `|` character between each
allele?**

``` text

Answer: Each of the SNVs have a '/'.
```

**Q6 - Check your notes or search what the difference between `/` or `|`
means in a VCF genotype. What does this mean?**


    Answer: The dash indicates that the genetic data is unphased; in other words, we do not know which allele comes from which parent. The pipe indicates that the alleles are listed in a specific chromosomal order; the former indicates the maternal and the latter indicates the paternal.

The reason HG002 has a \| operator in its genotypes is that it was
trio-binned, meaning we used sequences of HG003 and HG004 to determine
whether each allele falls in the paternal or maternal haplotype.

**Q7 - Write out HG002 haplotype 1 (paternal) and 2 (maternal) as a
sequence of alleles at SNV chr8:20,098,302, chr8:20,098,409, and
chr8:20,098,745.**

``` text

Haplotype 1/paternal:G, T, A

Haplotype 2/maternal: G, T, A
```

**Q8 - Can we determine the haplotypes of the other patients? Why or why
not?**


    Answer: We can, but on the condition that they have a '|' rather than a slash. So, hg001 and hg002 can be identified whereas all the others cannot.

**Q9 - Calculate the nucleotide diversity of chr8:20,098,302 in all
patients, not including the reference genome. Show your work.**


    Answer:
    Pi = 1/(n/2) Σ d
    diff pairs/all pairs
    24 / 91 -> 0.2637...

    14 haplotypes
    2 similarities for A and G
    G|A
    G|G
    G/G
    A/G
    G/G
    G/G
    G/G

### Step 2 - Analyzing the VCF in R

Run the code chunk below to load vcfR, an R package for VCF operations,
and import your VCF

``` r
library(vcfR)
```

    ## Warning: package 'vcfR' was built under R version 4.6.1

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.16.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 4.6.1

    ## Warning: package 'ggplot2' was built under R version 4.6.1

    ## Warning: package 'tibble' was built under R version 4.6.1

    ## Warning: package 'tidyr' was built under R version 4.6.1

    ## Warning: package 'readr' was built under R version 4.6.1

    ## Warning: package 'purrr' was built under R version 4.6.1

    ## Warning: package 'dplyr' was built under R version 4.6.1

    ## Warning: package 'stringr' was built under R version 4.6.1

    ## Warning: package 'forcats' was built under R version 4.6.1

    ## Warning: package 'lubridate' was built under R version 4.6.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.2.1     ✔ readr     2.2.0
    ## ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
    ## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
    ## ✔ purrr     1.2.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# import the vcf
vcf <- read.vcfR("hg0001_7_merged.vcf.gz")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 67
    ##   header_line: 68
    ##   variant count: 371213
    ##   column count: 16
    ## Meta line 67 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 371213
    ##   Character matrix gt cols: 16
    ##   skip: 0
    ##   nrows: 371213
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant 15000Processed variant 16000Processed variant 17000Processed variant 18000Processed variant 19000Processed variant 20000Processed variant 21000Processed variant 22000Processed variant 23000Processed variant 24000Processed variant 25000Processed variant 26000Processed variant 27000Processed variant 28000Processed variant 29000Processed variant 30000Processed variant 31000Processed variant 32000Processed variant 33000Processed variant 34000Processed variant 35000Processed variant 36000Processed variant 37000Processed variant 38000Processed variant 39000Processed variant 40000Processed variant 41000Processed variant 42000Processed variant 43000Processed variant 44000Processed variant 45000Processed variant 46000Processed variant 47000Processed variant 48000Processed variant 49000Processed variant 50000Processed variant 51000Processed variant 52000Processed variant 53000Processed variant 54000Processed variant 55000Processed variant 56000Processed variant 57000Processed variant 58000Processed variant 59000Processed variant 60000Processed variant 61000Processed variant 62000Processed variant 63000Processed variant 64000Processed variant 65000Processed variant 66000Processed variant 67000Processed variant 68000Processed variant 69000Processed variant 70000Processed variant 71000Processed variant 72000Processed variant 73000Processed variant 74000Processed variant 75000Processed variant 76000Processed variant 77000Processed variant 78000Processed variant 79000Processed variant 80000Processed variant 81000Processed variant 82000Processed variant 83000Processed variant 84000Processed variant 85000Processed variant 86000Processed variant 87000Processed variant 88000Processed variant 89000Processed variant 90000Processed variant 91000Processed variant 92000Processed variant 93000Processed variant 94000Processed variant 95000Processed variant 96000Processed variant 97000Processed variant 98000Processed variant 99000Processed variant 100000Processed variant 101000Processed variant 102000Processed variant 103000Processed variant 104000Processed variant 105000Processed variant 106000Processed variant 107000Processed variant 108000Processed variant 109000Processed variant 110000Processed variant 111000Processed variant 112000Processed variant 113000Processed variant 114000Processed variant 115000Processed variant 116000Processed variant 117000Processed variant 118000Processed variant 119000Processed variant 120000Processed variant 121000Processed variant 122000Processed variant 123000Processed variant 124000Processed variant 125000Processed variant 126000Processed variant 127000Processed variant 128000Processed variant 129000Processed variant 130000Processed variant 131000Processed variant 132000Processed variant 133000Processed variant 134000Processed variant 135000Processed variant 136000Processed variant 137000Processed variant 138000Processed variant 139000Processed variant 140000Processed variant 141000Processed variant 142000Processed variant 143000Processed variant 144000Processed variant 145000Processed variant 146000Processed variant 147000Processed variant 148000Processed variant 149000Processed variant 150000Processed variant 151000Processed variant 152000Processed variant 153000Processed variant 154000Processed variant 155000Processed variant 156000Processed variant 157000Processed variant 158000Processed variant 159000Processed variant 160000Processed variant 161000Processed variant 162000Processed variant 163000Processed variant 164000Processed variant 165000Processed variant 166000Processed variant 167000Processed variant 168000Processed variant 169000Processed variant 170000Processed variant 171000Processed variant 172000Processed variant 173000Processed variant 174000Processed variant 175000Processed variant 176000Processed variant 177000Processed variant 178000Processed variant 179000Processed variant 180000Processed variant 181000Processed variant 182000Processed variant 183000Processed variant 184000Processed variant 185000Processed variant 186000Processed variant 187000Processed variant 188000Processed variant 189000Processed variant 190000Processed variant 191000Processed variant 192000Processed variant 193000Processed variant 194000Processed variant 195000Processed variant 196000Processed variant 197000Processed variant 198000Processed variant 199000Processed variant 200000Processed variant 201000Processed variant 202000Processed variant 203000Processed variant 204000Processed variant 205000Processed variant 206000Processed variant 207000Processed variant 208000Processed variant 209000Processed variant 210000Processed variant 211000Processed variant 212000Processed variant 213000Processed variant 214000Processed variant 215000Processed variant 216000Processed variant 217000Processed variant 218000Processed variant 219000Processed variant 220000Processed variant 221000Processed variant 222000Processed variant 223000Processed variant 224000Processed variant 225000Processed variant 226000Processed variant 227000Processed variant 228000Processed variant 229000Processed variant 230000Processed variant 231000Processed variant 232000Processed variant 233000Processed variant 234000Processed variant 235000Processed variant 236000Processed variant 237000Processed variant 238000Processed variant 239000Processed variant 240000Processed variant 241000Processed variant 242000Processed variant 243000Processed variant 244000Processed variant 245000Processed variant 246000Processed variant 247000Processed variant 248000Processed variant 249000Processed variant 250000Processed variant 251000Processed variant 252000Processed variant 253000Processed variant 254000Processed variant 255000Processed variant 256000Processed variant 257000Processed variant 258000Processed variant 259000Processed variant 260000Processed variant 261000Processed variant 262000Processed variant 263000Processed variant 264000Processed variant 265000Processed variant 266000Processed variant 267000Processed variant 268000Processed variant 269000Processed variant 270000Processed variant 271000Processed variant 272000Processed variant 273000Processed variant 274000Processed variant 275000Processed variant 276000Processed variant 277000Processed variant 278000Processed variant 279000Processed variant 280000Processed variant 281000Processed variant 282000Processed variant 283000Processed variant 284000Processed variant 285000Processed variant 286000Processed variant 287000Processed variant 288000Processed variant 289000Processed variant 290000Processed variant 291000Processed variant 292000Processed variant 293000Processed variant 294000Processed variant 295000Processed variant 296000Processed variant 297000Processed variant 298000Processed variant 299000Processed variant 300000Processed variant 301000Processed variant 302000Processed variant 303000Processed variant 304000Processed variant 305000Processed variant 306000Processed variant 307000Processed variant 308000Processed variant 309000Processed variant 310000Processed variant 311000Processed variant 312000Processed variant 313000Processed variant 314000Processed variant 315000Processed variant 316000Processed variant 317000Processed variant 318000Processed variant 319000Processed variant 320000Processed variant 321000Processed variant 322000Processed variant 323000Processed variant 324000Processed variant 325000Processed variant 326000Processed variant 327000Processed variant 328000Processed variant 329000Processed variant 330000Processed variant 331000Processed variant 332000Processed variant 333000Processed variant 334000Processed variant 335000Processed variant 336000Processed variant 337000Processed variant 338000Processed variant 339000Processed variant 340000Processed variant 341000Processed variant 342000Processed variant 343000Processed variant 344000Processed variant 345000Processed variant 346000Processed variant 347000Processed variant 348000Processed variant 349000Processed variant 350000Processed variant 351000Processed variant 352000Processed variant 353000Processed variant 354000Processed variant 355000Processed variant 356000Processed variant 357000Processed variant 358000Processed variant 359000Processed variant 360000Processed variant 361000Processed variant 362000Processed variant 363000Processed variant 364000Processed variant 365000Processed variant 366000Processed variant 367000Processed variant 368000Processed variant 369000Processed variant 370000Processed variant 371000Processed variant: 371213
    ## All variants processed

``` r
# peak at the vcf
head(vcf)
```

    ## [1] "***** Object of class 'vcfR' *****"
    ## [1] "***** Meta section *****"
    ## [1] "##fileformat=VCFv4.2"
    ## [1] "##FILTER=<ID=PASS,Description=\"All filters passed\">"
    ## [1] "##fileDate=20160824"
    ## [1] "##CL=vcffilter -i filtered-phase-transfer.vcf.gz -o - --javascript \" [Truncated]"
    ## [1] "##RUN-ID=16dacf15-fdc9-4199-84bd-723ea8bcddef"
    ## [1] "##contig=<ID=chr1,length=248956422>"
    ## [1] "First 6 rows."
    ## [1] 
    ## [1] "***** Fixed section *****"
    ##      CHROM  POS      ID REF ALT QUAL FILTER
    ## [1,] "chr8" "221230" NA "T" "G" "50" "PASS"
    ## [2,] "chr8" "223423" NA "T" "C" "50" "PASS"
    ## [3,] "chr8" "225124" NA "G" "C" "50" "PASS"
    ## [4,] "chr8" "226273" NA "G" "T" "50" "PASS"
    ## [5,] "chr8" "226818" NA "G" "A" "50" "PASS"
    ## [6,] "chr8" "239493" NA "G" "A" "50" "PASS"
    ## [1] 
    ## [1] "***** Genotype section *****"
    ##      FORMAT                         HG001              
    ## [1,] "GT:PS:DP:ADALL:AD:GQ"         "./.:.:.:.:.:."    
    ## [2,] "GT:DP:ADALL:AD:GQ:IGT:IPS:PS" "./.:.:.:.:.:.:.:."
    ## [3,] "GT:PS:DP:ADALL:AD:GQ"         "./.:.:.:.:.:."    
    ## [4,] "GT:PS:DP:ADALL:AD:GQ"         "./.:.:.:.:.:."    
    ## [5,] "GT:DP:ADALL:AD:GQ:IGT:IPS:PS" "./.:.:.:.:.:.:.:."
    ## [6,] "GT:DP:ADALL:AD:GQ:IGT:IPS:PS" "./.:.:.:.:.:.:.:."
    ##      HG002                              HG003               HG004              
    ## [1,] "./.:.:.:.:.:."                    "./.:.:.:.:.:."     "./.:.:.:.:.:."    
    ## [2,] "0/1:900:259,241:0,0:594:0/1:.:."  "./.:.:.:.:.:.:.:." "./.:.:.:.:.:.:.:."
    ## [3,] "./.:.:.:.:.:."                    "./.:.:.:.:.:."     "./.:.:.:.:.:."    
    ## [4,] "./.:.:.:.:.:."                    "./.:.:.:.:.:."     "./.:.:.:.:.:."    
    ## [5,] "0/1:1056:280,285:0,0:495:0/1:.:." "./.:.:.:.:.:.:.:." "./.:.:.:.:.:.:.:."
    ## [6,] "0/1:1036:272,293:0,0:495:0/1:.:." "./.:.:.:.:.:.:.:." "./.:.:.:.:.:.:.:."
    ##      HG005              
    ## [1,] "./.:.:.:.:.:."    
    ## [2,] "./.:.:.:.:.:.:.:."
    ## [3,] "./.:.:.:.:.:."    
    ## [4,] "./.:.:.:.:.:."    
    ## [5,] "./.:.:.:.:.:.:.:."
    ## [6,] "./.:.:.:.:.:.:.:."
    ## [1] "First 6 columns only."
    ## [1] 
    ## [1] "Unique GT formats:"
    ## [1] "GT:PS:DP:ADALL:AD:GQ"         "GT:DP:ADALL:AD:GQ:IGT:IPS:PS"
    ## [1]

``` r
# create a genotype matrix from vcf
gt <- extract.gt(vcf, element = "GT")

# peak at the matrix
head(gt)
```

    ##               HG001 HG002 HG003 HG004 HG005 INTEGRATION 7:INTEGRATION
    ## chr8_221230_1 NA    NA    NA    NA    NA    NA          "0/1"        
    ## chr8_223423_2 NA    "0/1" NA    NA    NA    NA          NA           
    ## chr8_225124_3 NA    NA    NA    NA    NA    "1/1"       NA           
    ## chr8_226273_4 NA    NA    NA    NA    NA    NA          "0/1"        
    ## chr8_226818_5 NA    "0/1" NA    NA    NA    NA          NA           
    ## chr8_239493_6 NA    "0/1" NA    NA    NA    NA          NA

``` r
colnames(gt)
```

    ## [1] "HG001"         "HG002"         "HG003"         "HG004"        
    ## [5] "HG005"         "INTEGRATION"   "7:INTEGRATION"

1.  Calculate heterozygosity:

Finish the code below under the commented area for calculating
heterozygosity per individual This code chunk starts by counting the
number of heterozygous genotypes in the VCF per individual. The
proportion of heterozygosity across the genome of each individual is
then calculated as…

$$ \text{Het} = \frac{\text{Het calls}}{\text{callable genotypes}} $$
Callable genotypes are those with a genotype quality score greater than
some threshold, like \> 30.

A vector with the number of callable genotypes per individual is created
below. Use the variables `het$sites` and `callable` to calculate
heterozygosity as the variable `het$prop_het`.

``` r
# calculate heterozygosity for each individual
het_sites <- colSums(gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0", na.rm = TRUE)

# print results
het <- as.data.frame(het_sites)

# this is a vector of the number of genotype-callable sites
callable <-     c(126578828, # hg001 
                  130612178, 123427582, 123931918, # hg002-4 trio
                  125738207, 129846161, 130150817) # hg005-7 trio

#######################
### YOUR CODE BELOW ###
#######################

# Use the variables `het$sites` and `callable` to calculate heterozygosity as the variable `het$prop_het`.
het$prop_het <- (het_sites)/(callable)


# print results
het
```

    ##               het_sites     prop_het
    ## HG001            110824 0.0008755335
    ## HG002            124370 0.0009522083
    ## HG003            112111 0.0009083140
    ## HG004            107232 0.0008652493
    ## HG005             98609 0.0007842405
    ## INTEGRATION      100000 0.0007701421
    ## 7:INTEGRATION    106850 0.0008209706

2.  Perform a principal components analysis of genetic variation across
    all individuals.

Run the code chunk below to convert genotypes to minor allele
frequencies: the number of alleles of the minor allele per genotype.

``` r
# convert genotypes to minor allele frequencies
gt[gt == "0/0" | gt == "0|0"] <- 0
gt[gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0"] <- 1
gt[gt == "1/1" | gt == "1|1"] <- 2

# convert to numeric matrix
gt <- apply(gt, 2, as.numeric)
```

    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion
    ## Warning in apply(gt, 2, as.numeric): NAs introduced by coercion

``` r
# transpose so rows = samples and columns = SNPs
gt <- t(gt)

# remove SNPs with missing genotypes
gt <- gt[, colSums(is.na(gt)) == 0]
```

Finish the code below to perform a PCA using the function `prcomp()`
input with the genotype data `gt`. Withint the function set
`center = TRUE`, and `scale. = FALSE`.

Review the prcomp() manual here to determine what the center and scale.
parameters do:
<https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/prcomp>

``` r
#######################
### YOUR CODE BELOW ###
#######################

# run PCA using the function prcomp()
# input the data frame gt to the function
# set center = TRUE
# set scale. = FALSE

pca <- prcomp(gt, center= TRUE, scale.=FALSE)
```

Run the code chunk below to calculate how much genetic variation is
explained by principal components 1, 2 and 3. Then, the code puts that
information into a data frame.

``` r
# calculate variance explained
variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# create plotting data frame
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3]
)
```

Finish the code below to plot your PCA analysis. First, finish the
incomplete ggplot code to plot samples on PC1 and PC2. Then, copy your
completed ggplot code, paste it, and edit it to plot PC1 and PC3.

``` r
#############################
### FINISH THE CODE BELOW ###
#############################

# plot PC1 and PC2
ggplot(pca_df, aes(x=PC1, y=PC2, label = Sample))+
  
  geom_point() +
  geom_text(vjust = -0.8, size = 3) +
  labs(
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC2 (", round(variance[2], 2), "%)")) +
  theme_classic()
```

![](bi252_wk9_lab_notebook_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# plot PC1 and PC3 by using your by copy/pasting your code above but adjusting your Y axis to PC3 and adjusting your axis labels in labs()
ggplot(pca_df, aes(x=PC1, y=PC3, label = Sample))+
  
  geom_point() +
  geom_text(vjust = -0.8, size = 3) +
  labs(
    x = paste0("PC1 (", round(variance[1], 2), "%)"),
    y = paste0("PC3 (", round(variance[3], 2), "%)")) +
  theme_classic()
```

![](bi252_wk9_lab_notebook_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

**Q10 - What principal component separated individuals based on
ancestry?**

``` text

Answer: PC1 separated individuals based on ancestry, as can be seen in the graphs above.  
```

**Q11 - If you consider the proportions of variation explained by
components 1 - 3, was a majority of genetic variation explained by known
differences in ancestry between individuals?**

``` text

Answer: For the most part, a majority of the genetic variation can be explained by known differences in individuals. Notably, though, the HG002/3/4 family has more variation than the HG005/6/7 family.
```

### Step 3 - Calculate number of breeding individuals (Ne) in each patient’s ancestral population

Recall that Ne can be calculated from heterozygosity (He) with:

$$N_e = \frac{H_e}{4\mu}$$

The human mutation rate is 1.25e-8 (or 1.25 x 10^-8). With this
information, finish the code below to report Ne of each individual’s
ancestral populations.

``` r
# using het$prop_het, add a new variable to the het df that calculates Ne (effective population size)

#############################
### FINISH THE CODE BELOW ###
#############################

het$Ne <-  het$prop_het/(4*(1.25 * 10^-8))

# report Ne
het$Ne
```

    ## [1] 17510.67 19044.17 18166.28 17304.99 15684.81 15402.84 16419.41

**Q12 - Ne is interpreted as the average number of breeding individuals
within a population per generation over the ancestral history of an
individual or population. Are the values that you’re seeing higher or
lower than you would have expected? Why?**

``` text

Answer: The values are different than expected. For example, it was expected that the number of breeding indiviudals in the Han Chinese groups would be larger than the Ashkenazi Jewish group. Now looking back, though, it makes more sense. The three individuals of the Han Chinese ancestry had closer relation to eachother than the Ashkenazi Jewish family.
```

**Q13 - Two populations of the same census size and different
heterozygosities will have different Ne. Humans have relatively low
heterozygosity and high inbreeding compared to most animals. Does this
information change your expectations about the Ne of human populations?
Why?**

``` text

Answer: Not exactly. We do not have a good reference point for what an Ne value would look like--controlled for population size--when there is high versus low inbreeding. Obviously, the former would have a smaller Ne value and the latter would have a higher one, but to what degree the Ne would be affected we do not know.
```

### Step 4 - Test the hypothesis that migratory distance from ancestral Homo sapiens population explains genetic diversity

Finish the code below to add km of distance from Addis Ababa, Ethiopia,
as a variable to the `het` data frame called `het$aa_dist`

Then, write ggplot code to plot aa_dist on the x axis for each patient
and prop_het on the y axis. Within ggplot(), set data = het.

``` r
# add distance from addis ababa (km) as a variable to the het data frame. reference the distances at the top of this notebook.

het$aa_dist <- c(5600, 2600, 2600, 2600, 8300, 8300, 8300)


# make a figure with ggplot and geom_point() with addis ababa distance on the x axis and prop_het on the y axis. use the het data frame as the input.

#######################
### YOUR CODE BELOW ###
#######################

ggplot(het, aes(x=aa_dist, y=prop_het, label = rownames(pca$x)))+
  geom_point() +
  geom_text(vjust = -0.8, size = 3) +
  labs(
    x = paste0("Distance"),
    y = paste0("Prop_het")) +
  theme_classic()
```

![](bi252_wk9_lab_notebook_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

**Q14 - Write a hypothesis that this analysis would test using an if,
then sentence**

``` text

Hypothesis: If distance in km from Addis Ababa is higher, then Ne will be lower.
```

**Q15 - Do your results support or reject your hypothesis?**

``` text

Answer: Yes, there appears to be a correlation as we described. This is likely due to the founder effect, which stipulates that as smaller parts of a population form new populations (through geographical movement, for example), the genetic diversity becomes smaller. In the new population, there is a smaller number of individuals in the breeding population which is ultimately what causes the lowered genetic diversity. In the case of our specific data, this explains the genetic closeness of the Han Chinese group compared to the patients of other ancestry. The people occupying the Han Chinese region had to go through more founder effect events.
```

## Submission

Knit this .Rmd file. Upload the knitted .md file and associated figures
to your lab Github repository.
