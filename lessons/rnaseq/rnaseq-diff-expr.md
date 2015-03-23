---
layout: page
---



# RNA-seq: differential gene expression analysis

This is an introduction to RNAseq analysis involving reading in count data from an RNAseq experiment, exploring the data using base R functions and then analysis with the DESeq2 package.

## Install and load packages

First, we'll need to install some add-on packages. Most generic R packages are hosted on the Comprehensive R Archive Network (CRAN, <http://cran.us.r-project.org/>). To install one of these packages, you would use `install.packages("packagename")`. You only need to install a package once, then load it each time using `library(packagename)`. Let's install the **gplots** and **calibrate** packages.


```r
install.packages("gplots")
install.packages("calibrate")
```

Bioconductor packages work a bit differently, and are not hosted on CRAN. Go to <http://bioconductor.org/> to learn more about the Bioconductor project. To use any Bioconductor package, you'll need a few "core" Bioconductor packages. Run the following commands to (1) download the installer script, and (2) install some core Bioconductor packages. You'll need internet connectivity to do this, and it'll take a few minutes, but it only needs to be done once.


```r
# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages.
# This requires internet connectivity and will take some time!
biocLite()
```

To install specific packages, first download the installer script if you haven't done so, and use `biocLite("packagename")`. This only needs to be done once then you can load the package like any other package. Let's download the [DESeq2 package](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html):


```r
# Do only once
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

Now let's load the packages we'll use:


```r
library(DESeq2)
library(gplots)
library(calibrate)
```

Bioconductor packages usually have great documentation in the form of *vignettes*. For a great example, take a look at the [DESeq2 vignette for analyzing count data](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).

## Introduction and data import

Analyzing an RNAseq experiment begins with sequencing reads. These are aligned to a reference genome, then the number of reads mapped to each gene can be counted.

The data for this tutorial comes from a PLOS ONE paper, [Genome-Wide Transcriptional Profiling of Skin and Dorsal Root Ganglia after Ultraviolet-B-Induced Inflammation](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0093338), and the raw data can be downloaded from [Gene Expression Omnibus database (GEO)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54413).

This data has already been downloaded and aligned to the human genome. The command line tool [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) was used to count reads mapped to human genes from the [Ensembl annotation](http://www.ensembl.org/info/data/ftp/index.html).

The output from this tool is provided in the `counts.txt` file in the `data` directory. Have a look at this file in the shell, using `head`.

First, set your working directory to the top level of the RNA-seq course. Import the data into R as a `data.frame` and examine it again. You can set the arguments of `read.table` to import the first row as a header giving the column names, and the first column as row names.


```r
# Read in the data
mycounts <- read.table("data/counts.txt")
head(mycounts)
mycounts <- read.table("data/counts.txt", header = TRUE)
head(mycounts)
mycounts <- read.table("data/counts.txt", header = TRUE, row.names = 1)
head(mycounts)
rownames(mycounts)
colnames(mycounts)
class(mycounts)
```

The data.frame contains information about genes (one gene per row) with the gene positions in the first five columns and then information about the number of reads aligning to the gene in each experimental sample. There are three replicates for control (column names starting with "ctl") and three for samples treated with ultraviolet-B light (starting "uvb"). We don't need the information on gene position for this analysis, just the counts for each gene and sample, so we can remove it from the data frame.


```r
# First, take a look
options(max.print = 80)
head(mycounts)
colnames(mycounts)

# One way to do it: by indexes
mycounts[, -(1:5)]

# Another way to do it, by arg to subset
subset(mycounts, select = c(-Chr, -Start, -End, -Strand, -Length))
subset(mycounts, select = -c(Chr, Start, End, Strand, Length))

# Reassign
mycounts <- subset(mycounts, select = -c(Chr, Start, End, Strand, Length))

# Take another look
head(mycounts)
colnames(mycounts)
```

---

**EXERCISE 1**

There's an R function called `rowSums()` that calculates the sum of each row in a numeric matrix, like the count matrix we have here, and it returns a vector. There's also a function called `which.max()` that determines the index of the maximum value in a vector. Here's an example of it in action


```r
# Create a temporary data frame called fake
fake <- data.frame(row.names = c("GeneA", "GeneB", "GeneC"),
                  samp1=c(20,50,40), samp2=c(30,70,50))

# This is what it looks like
fake
```

```
##       samp1 samp2
## GeneA    20    30
## GeneB    50    70
## GeneC    40    50
```

```r
# Get the rowSums
rowSums(fake)
```

```
## GeneA GeneB GeneC 
##    50   120    90
```

```r
# Get the index of the maximum total value
which.max(rowSums(fake))
```

```
## GeneB 
##     2
```

```r
# Store that index
topGene <- which.max(rowSums(fake))

# Get that row, and all the columns
fake[topGene, ]
```

```
##       samp1 samp2
## GeneB    50    70
```


0. Find the gene with the highest expression across all samples -- remember, each row is a gene.
0. Extract the expression data for this gene for all samples.
0. In which sample does it have the highest expression?
0. What is the function of the gene? Can you suggest why this is the top expressed gene?



---

## Data investigation using base R

We can investigate this data a bit more using some of the basic R functions before going on to use more sophisticated analysis tools.

First make a copy of the data, because we'll need it later. We will work on the copy. We will calculate the mean for each gene for each condition and plot them.


```r
mc2 <- mycounts  #make a copy

# get Control columns
colnames(mc2)

# grep searches for matches to a pattern. Get help with ?grep

# get the indexes for the controls
grep("ctl", colnames(mc2))
ctlCols <- grep("ctl", colnames(mc2))
mc2[, ctlCols]

# use the rowMeans function
mc2$ctlmean <- rowMeans(mc2[, ctlCols])
head(mc2)

# same for uvb
uvbCols <- grep("uvb", colnames(mc2))
mc2$uvbmean <- rowMeans(mc2[, uvbCols])

head(mc2)
```

---

**EXERCISE 2**

Using the `subset()` function, print out all the columns where the control mean does not equal 0 **and** where the UVB mean does not equal zero.

Bonus: code golf -- use the fewest characters to get the same solution.



---

**EXERCISE 3**

1. Plot the mean expression of each gene in control against the UVB sample mean. Are there any outliers?
2. How could you make this plot more informative and look more professional? Hint: try plotting on the log scale and using a different point character.


```r
plot(mc2$ctlmean, mc2$uvbmean)
with(mc2, plot(log10(ctlmean), log10(uvbmean), pch = 16))
with(mc2, plot(ctlmean, uvbmean, log = "xy", pch = 16))
```

---

## Poor man's differential gene expression

We can find candidate differentially expressed genes by looking for genes with a large change between control and UVB samples. A common threshold used is log2 fold change more than 2 or less than -2. We will calculate log2 fold change for all the genes and colour the genes with log2 fold change of more than 2 or less than -2 on the plot.

First, check for genes with a mean expression of 0. Putting zeroes into the log2 fold change calculation will produce NAs, so we might want to remove these genes. Note: this is for mathematical reasons, although different software may produce different results when you try to do `log2(0)`.

In R, `TRUE` and `FALSE` can be represented as `1` and `0`, respectively. When we call `sum(mc2$ctlmean > 0)`, we're really asking, "how many genes have a mean above 0 in the control group?"


```r
sum(mc2$ctlmean > 0)
sum(mc2$uvbmean > 0)
```

Now, let's subset the data and keep things where either the control or UVB group means are greater than zero.


```r
nrow(mc2)
mc2 <- subset(mc2, mc2$ctlmean > 0 | mc2$uvbmean > 0)
nrow(mc2)
head(mc2)
```

Mathematically things work out better for us when we test things on the log scale. On the absolute scale, upregulation goes from 1 to infinity, while downregulation is bounded by 0 and 1. On the log scale, upregulation goes from 0 to infinity, and downregulation goes from 0 to negative infinity. Let's compute a log-base-2 of the fold change.

When we do this we'll see some `Inf` and `-Inf` values. This is what happens when we take `log2(Inf)` or `log2(0)`.


```r
# calculate the log2 fold change
mc2$log2FC <- log2(mc2$uvbmean/mc2$ctlmean)
head(mc2)

# see how many are up, down, or both (using the absolute value function)
sum(mc2$log2FC > 2)
sum(mc2$log2FC < -2)
sum(abs(mc2$log2FC) > 2)
```

We can few just the "outliers" with `subset(mc2, abs(log2FC)>2)`, which gets us just the things that have a large fold change in either direction. Let's plot these.


```r
with(mc2, plot(ctlmean, uvbmean, log = "xy", pch = 16))

# Use the points function to add to the same plot
mc2de <- subset(mc2, abs(log2FC) > 2)
with(mc2de, points(ctlmean, uvbmean, pch = 16, col = "red"))
```

What do you notice about the positions of the outliers on these plots? How would you interpret this? What are some of the problems with this simple approach?

## DESeq2 analysis

DESeq2 is an R package for analyzing count-based NGS data like RNA-seq. It is available from [Bioconductor](http://www.bioconductor.org/). Bioconductor is a project to provide tools for analysing high-throughput genomic data including RNA-seq, ChIP-seq and arrays. You can explore Bioconductor packages [here](http://www.bioconductor.org/packages/release/BiocViews.html#___Software).

Just like R packages from CRAN, you only need to install Bioconductor packages once, then load them every time you start a new R session.


```r
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```


```r
library("DESeq2")
citation("DESeq2")
```

It requires the count data to be in matrix form, and an additional dataframe describing sample metadata. Notice that the **colnames of the countdata** match the **rownames of the metadata*. If you read it in without specifying which column contains the row names, the row names will be blank. We have to specify column 1 contains the row names.


```r
mycoldata <- read.csv("data/coldata.csv", header = TRUE)
mycoldata
mycoldata <- read.csv("data/coldata.csv", header = TRUE, row.names = 1)
mycoldata
```

DESeq works on a particular type of object called a DESeqDataSet. The DESeqDataSet is a single object that contains input values, intermediate calculations like how things are normalized, and all results of a differential expression analysis. You can construct a DESeqDataSet from a count matrix, a metadata file, and a formula indicating the design of the experiment.


```r
dds <- DESeqDataSetFromMatrix(countData = mycounts, colData = mycoldata, design = ~condition)
dds
```

Next, let's run the DESeq pipeline on the dataset, and reassign the resulting object back to the same variable. Before we start, `dds` is a bare-bones DESeqDataSet. The `DESeq()` function takes a DESeqDataSet and returns a DESeqDataSet, but with lots of other information filled in (normalization, results, etc). Here, we're running the DESeq pipeline on the `dds` object, and reassigning the whole thing back to `dds`, which will now be a DESeqDataSet populated with results.


```r
dds <- DESeq(dds)
```

Now, let's use the `results()` function to pull out the results from the `dds` object. Let's re-order by the adjusted p-value.


```r
# Get differential expression results
res <- results(dds)
res

# Order by adjusted p-value
order(res$padj)
res[order(res$padj), ]
res <- res[order(res$padj), ]
res
```

Combine DEseq results with the original counts data. Write significant results to a file.


```r
sig <- subset(res, padj < 0.05)
dir.create("results")
write.csv(sig, file = "results/sig.csv")  # tab delim data
```

You can open this file in Excel or any text editor (try it now).

## Data Visualization

We can also do some exploratory plotting of the data.

The differential expression analysis above operates on the raw (normalized) count data. But for visualizing or clustering data as you would with a microarray experiment, you ned to work with transformed versions of the data. First, use a *regularlized log* transofmration while re-estimating the dispersion ignoring any information you have about the samples (`blind=TRUE`). Perform a principal components analysis and hierarchical clustering.


```r
# Transform
rld <- rlogTransformation(dds)

# Principal components analysis
plotPCA(rld, intgroup = "condition")

# Hierarchical clustering analysis let's get the actual values for the first
# few genes
head(assay(rld))
## now transpose those
t(head(assay(rld)))
## now get the sample distances from the transpose of the whole thing
dist(t(assay(rld)))
sampledist <- dist(t(assay(rld)))
plot(hclust(sampledist))
```

Let's plot a heatmap.


```r
# ?heatmap for help
sampledist
as.matrix(sampledist)
sampledistmat <- as.matrix(sampledist)
heatmap(sampledistmat)
```

That's a horribly ugly default. You can change the built-in heatmap function, but others are better.


```r
# better heatmap with gplots
library("gplots")
heatmap.2(sampledistmat)
heatmap.2(sampledistmat, key = FALSE, trace = "none")
colorpanel(10, "black", "white")
heatmap.2(sampledistmat, col = colorpanel(64, "black", "white"), key = FALSE, 
    trace = "none")
heatmap.2(sampledistmat, col = colorpanel(64, "steelblue", "white"), key = FALSE, 
    trace = "none")
heatmap.2(sampledistmat, col = colorpanel(64, "red", "white", "blue"), key = FALSE, 
    trace = "none")
```

What about a histogram of the p-values?


```r
# Examine plot of p-values
hist(res$pvalue, breaks = 50, col = "grey")
```

Let's plot an MA-plot. This shows the fold change versus the overall expression values.


```r
with(res, plot(baseMean, log2FoldChange, pch = 16, cex = 0.5, log = "x"))
with(subset(res, padj < 0.05), points(baseMean, log2FoldChange, col = "red", 
    pch = 16))

# optional: label the points with the calibrate package. see ?textxy for
# help
library("calibrate")
res$Gene <- rownames(res)
head(res)
with(subset(res, padj < 0.05), textxy(baseMean, log2FoldChange, labs = Gene, 
    cex = 1, col = "red"))
```

Let's create a volcano plot.


```r
par(pch = 16)
with(res, plot(log2FoldChange, -log10(pvalue), main = "Volcano plot"))
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), col = "red"))
with(subset(res, abs(log2FoldChange) > 2), points(log2FoldChange, -log10(pvalue), 
    col = "orange"))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), points(log2FoldChange, 
    -log10(pvalue), col = "green"))
# Add legend
legend("topleft", legend = c("FDR<0.05", "|LFC|>2", "both"), pch = 16, col = c("red", 
    "orange", "green"))
# Label points
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), textxy(log2FoldChange, 
    -log10(pvalue), labs = Gene, cex = 1))
```

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
##  [1] calibrate_1.7.2           MASS_7.3-39              
##  [3] gplots_2.16.0             DESeq2_1.6.3             
##  [5] RcppArmadillo_0.4.650.1.1 Rcpp_0.11.5              
##  [7] GenomicRanges_1.18.4      GenomeInfoDb_1.2.4       
##  [9] IRanges_2.0.1             S4Vectors_0.4.0          
## [11] BiocGenerics_0.12.1       knitr_1.9                
## [13] BiocInstaller_1.16.1     
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3      annotate_1.44.0      AnnotationDbi_1.28.1
##  [4] base64enc_0.1-2      BatchJobs_1.5        BBmisc_1.9          
##  [7] Biobase_2.26.0       BiocParallel_1.0.2   bitops_1.0-6        
## [10] brew_1.0-6           caTools_1.17.1       checkmate_1.5.1     
## [13] cluster_2.0.1        codetools_0.2-10     colorspace_1.2-4    
## [16] DBI_0.3.1            digest_0.6.8         evaluate_0.5.5      
## [19] fail_1.2             foreach_1.4.2        foreign_0.8-63      
## [22] formatR_1.0          Formula_1.2-0        gdata_2.13.3        
## [25] genefilter_1.48.1    geneplotter_1.44.0   ggplot2_1.0.0       
## [28] grid_3.1.2           gtable_0.1.2         gtools_3.4.1        
## [31] Hmisc_3.15-0         iterators_1.0.7      KernSmooth_2.23-14  
## [34] labeling_0.3         lattice_0.20-30      latticeExtra_0.6-26 
## [37] locfit_1.5-9.1       munsell_0.4.2        nnet_7.3-9          
## [40] plyr_1.8.1           proto_0.3-10         RColorBrewer_1.1-2  
## [43] reshape2_1.4.1       rpart_4.1-9          RSQLite_1.0.0       
## [46] scales_0.2.4         sendmailR_1.2-1      splines_3.1.2       
## [49] stringr_0.6.2        survival_2.38-1      tools_3.1.2         
## [52] XML_3.98-1.1         xtable_1.7-4         XVector_0.6.0
```


## Going further

* After the course, download the [Integrative Genome Viewer](http://www.broadinstitute.org/igv/) from the Broad Institute. Download all your .bam files from your AWS instance, and load them into IGV. Try navigating to regions around differentially expressed genes to view how reads map to genes differently in the controls versus the irradiated samples.
* Can you see any genes where differential expression is likely attributable to a specific isoform?
* Do you see any instances of differential exon usage? You can investigate this formally with the [DEXSeq](http://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html) package.
* Read about pathway analysis with [GOSeq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html) or [SeqGSEA](http://www.bioconductor.org/packages/release/bioc/html/SeqGSEA.html) - tools for gene ontology analysis and gene set enrichment analysis using next-generation sequencing data.
* Read about multifactor designs in the [DESeq2 vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) for cases where you have multiple variables of interest (e.g. irradiated vs controls in multiple tissue types).

***After the course, make sure you stop any running AWS instances.***
