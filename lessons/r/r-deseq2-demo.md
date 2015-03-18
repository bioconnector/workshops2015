---
layout: page
---



# Demo: Using R and Bioconductor for RNA-seq data analysis

In this demo we'll analyze some publicly available RNA-seq gene expression data using R and bioinformatics-focused R packages in [Bioconductor](http://bioconductor.org/). This demo assumes some familiarity with R (functions, functions, vectors, creating variables, getting help, subsetting, data frames, plotting, and reading/writing files). The starting point for this analysis is a *count matrix* - RNA-seq data has already been cleaned, aligned, and counted to the gene level, and a *metadata file* with information about the samples. The `rownames` of the metadata must match the `colnames` of the count data.

We'll be using the DESeq2 package that you've already installed (see the setup instructions). 


```r
library(DESeq2)
```

Bioconductor packages usually have great documentation in the form of *vignettes*. For a great example, take a look at the [DESeq2 vignette for analyzing count data](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).




## Publicly available RNA-seq data

Now, let's analyze some publicly available gene expression (RNA-seq) data. NCBI Gene Expression Omnibus (<http://www.ncbi.nlm.nih.gov/geo/>) is an international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community. Many publishers require gene expression data be submitted to GEO and made publicly available before publication. You can learn a lot more about GEO by reading their [overview](http://www.ncbi.nlm.nih.gov/geo/info/overview.html) and [FAQ](http://www.ncbi.nlm.nih.gov/geo/info/faq.html) pages. At the time of this writing, GEO hosts over 45,000 studies comprising over 1,000,000 samples on over 10,000 different technology platforms.

In this demonstration, we're going to be using data from GEO Series accession number GSE18508. You can enter this number in the search box on the GEO homepage, or use [this direct link](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). In this study, the authors performed RNA-sequencing of mRNA from *Drosophila melanogaster* S2-DRSC cells that have been RNAi depleted of mRNAs encoding RNA binding proteins and splicing factors. This was done as part of the modENCODE project, published in: Brooks AN et al. Conservation of an RNA regulatory map between Drosophila and mammals. *Genome Res* 2011 Feb;21(2):193-202.

You can go to the [GEO accession page](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508) and download all the raw data at the bottom under supplementary data. However, this dataset is nearly 100GB. What I have in this repository is a spreadsheet containing a matrix of gene counts, where genes are in rows and samples are in columns. That is, this data has already been aligned and counted, and the number in the cell is the number of RNA-seq reads that mapped to that gene for that sample. The value in the *i*-th row and the *j*-th column of the matrix tells how many reads have been mapped to gene *i* in sample *j*. To do these steps yourself, you would need to align reads to the genome (e.g., using [STAR](https://code.google.com/p/rna-star/)) and count reads mapping to features (e.g., using a GTF file from Ensembl and a tool like [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)).

Note: much of this was adapted from the [DESeq2 package vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).

## Load the data

First, make sure your working directory (folder) is set to wherever you saved the data for this lesson. If you downloaded the [code repository for this lesson](https://github.com/bioconnector/workshops/archive/master.zip), that directory is `workshops/lessons/r`. This way, you can reference the data using a *relative path* (e.g. `data/pasilla_counts.csv`) instead of an *absolute path* (e.g. `C:/Users/name/downloads/workshops/lessons/data/pasilla_counts.csv`). You can check where you are with `getwd()`.


```r
getwd()
```

Let's look at the data we're about to load using Excel. The count data is stored in `data/pasilla_counts.csv` and the metadata is stored in `data/pasilla_metadata.csv`. 

Notice how cell A1 is empty. The first row corresponds to the `colnames`, and the first column corresponds to the `rownames`. When we read this data into R we need to be explicit about telling R that file has a column header (`header=TRUE`) and that the `rownames` for the data frame come from the first column (`row.names=1`).


```r
# Load the count data
pasillacounts <- read.csv("data/pasilla_counts.csv", header = TRUE, row.names = 1)
head(pasillacounts)

# Load the sample metadata
pasillameta <- read.csv("data/pasilla_metadata.csv", header = TRUE, row.names = 1)
pasillameta
```

Notice something here. The values that are printed across the top are the `colnames` and the values printed in the gutter to the side are the `rownames`. Notice how the **colnames of the countdata**, which are the names of our samples, match the **rownames of the metadata*. This is required to set up the data in a way for DESeq to work with it.

DESeq works on a particular type of object called a **DESeqDataSet**. The DESeqDataSet is a single object that contains input values, intermediate calculations like how things are normalized, and all results of a differential expression analysis. You can construct a DESeqDataSet from a count matrix, a metadata file, and a formula indicating the design of the experiment. The design formula expresses the variables which will be used in modeling. The formula should be a tilde (~) followed by the variables with plus signs between them. To do this, we run the `DESeqDataSetFromMatrix` function, giving it data frames containing the count matrix as well as the column (metadata) information, and the design formula.


```r
dds <- DESeqDataSetFromMatrix(countData = pasillacounts, colData = pasillameta, 
    design = ~condition)
dds
```

## Differential expression analysis

The standard differential expression analysis steps are wrapped into a single function, `DESeq`. This convenience function normalizes the library, estimates dispersion, and performs a statistical test under the negative binomial model.


```r
# Normalize, estimate dispersion, fit model
dds <- DESeq(dds)
```

The results are accessed using the function `results()`. The `mcols` function gives you more information about what the columns in the results tell you. 


```r
# Extract results
res <- results(dds)
head(res)
mcols(res)
```

Let's write out the results to a file. Then, let's create an MA-plot, showing the log2 fold change over the normalized counts for each gene, colored red if statistically significant (padj<0.1). Get more help with `?plotMA`.


```r
# Write out results
write.csv(res, file = "pasilla_results.csv")

# Create MA Plot
plotMA(dds)
```

## Data transformation and visualization

The differential expression analysis above operates on the raw (normalized) count data. But for visualizing or clustering data as you would with a microarray experiment, you ned to work with transformed versions of the data. First, use a *regularlized log* transformation then perform a principal components analysis.


```r
# Transform
rld <- rlogTransformation(dds)

# Principal components analysis
plotPCA(rld, intgroup = c("condition", "type"))
```

## Multifactor design

See the [DESeq2 vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) for more information on multifactor designs. This is how you would take library type into account (single-end vs paired-end sequencing).

## Record package and version info with `sessionInfo()`

The `sessionInfo()` prints version information about R and any attached packages. It's a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.


```r
sessionInfo()
```

```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] DESeq2_1.6.3              RcppArmadillo_0.4.650.1.1
##  [3] Rcpp_0.11.5               GenomicRanges_1.18.4     
##  [5] GenomeInfoDb_1.2.4        IRanges_2.0.1            
##  [7] S4Vectors_0.4.0           BiocGenerics_0.12.1      
##  [9] knitr_1.9                 BiocInstaller_1.16.1     
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3      annotate_1.44.0      AnnotationDbi_1.28.1
##  [4] base64enc_0.1-2      BatchJobs_1.5        BBmisc_1.9          
##  [7] Biobase_2.26.0       BiocParallel_1.0.2   brew_1.0-6          
## [10] checkmate_1.5.1      cluster_2.0.1        codetools_0.2-10    
## [13] colorspace_1.2-4     DBI_0.3.1            digest_0.6.8        
## [16] evaluate_0.5.5       fail_1.2             foreach_1.4.2       
## [19] foreign_0.8-63       formatR_1.0          Formula_1.2-0       
## [22] genefilter_1.48.1    geneplotter_1.44.0   ggplot2_1.0.0       
## [25] grid_3.1.2           gtable_0.1.2         Hmisc_3.15-0        
## [28] iterators_1.0.7      lattice_0.20-30      latticeExtra_0.6-26 
## [31] locfit_1.5-9.1       MASS_7.3-39          munsell_0.4.2       
## [34] nnet_7.3-9           plyr_1.8.1           proto_0.3-10        
## [37] RColorBrewer_1.1-2   reshape2_1.4.1       rpart_4.1-9         
## [40] RSQLite_1.0.0        scales_0.2.4         sendmailR_1.2-1     
## [43] splines_3.1.2        stringr_0.6.2        survival_2.38-1     
## [46] tools_3.1.2          XML_3.98-1.1         xtable_1.7-4        
## [49] XVector_0.6.0
```
