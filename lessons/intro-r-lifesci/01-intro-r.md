---
layout: page
---

# Introduction to R for Life Scientists



Link to slides: *Please check back after the workshop*

## Before coming

[Click here for installation instructions](../). You'll need to download R, RStudio, and install a few packages. The [README](README.md) will also give you some registration instructions.

## R basics

The first part of this workshop will demonstrate very basic functionality in r, including functions, functions, vectors, creating variables, getting help, subsetting, data frames, plotting, and reading/writing files.

### Basic operations

R can be used as a glorified calculator. Try typing this in directly into the console. Then start typing this into the editor, and save your script. Use the run button, or press `CMD`+`Enter` (`Ctrl`+`Enter` on Windows).


```r
2+2
5*4
2^3
```

R Knows order of operations.


```r
2+3*4/(5+3)*15/2^2+3*4^2
```

### Functions

R has built-in functions.


```r
# Notice that this is a comment.
# Anything behind a # is "commented out" and is not run.
sqrt(144)
log(1000)
```

Get help by typing a question mark in front of the function's name, or `help(functionname)`:

```
help(log)
?log
```

Note syntax highlighting when typing this into the editor. Also note how we pass *arguments* to functions. Finally, see how you can *next* one function inside of another (here taking the square root of the log-base-10 of 1000).


```r
log(1000)
log(1000, base=10)
sqrt(log(1000, base=10))
```

### Vectors

Let's create some numeric vectors. Vectors (aka "arrays" in Perl, "lists" in Python) are single *objects* containing an ordered collection of *elements*. A simple vector is a numeric vector, a single *object* containing several numbers. Here let's display a few vectors. We can also do vector arithmetic. When printing vectors to the screen that have lots of elements, notice that the bracketed number in the gutter of the output is just a counter indexing the number of elements in the vector.


```r
# Some simple numeric vectors:
1:5
6:10
1:5 + 6:10
1:100

# Get some help with the seq() function, then create a vector from 2 to 200 by 2s.
# Notice how the seq() function works -- the `to` argument will never be exceeded.
help(seq)
seq(from=2, to=200, by=4)
```

### Variables / objects

You can store values in a variable or object. Use the `<-` operator for assignment. `=` also will work, but `<-` is conventional and preferred. Objects should start with a letter and can include numbers and underscores. Named variables are objects containing whatever is assigned to them. Remember, *everything is an object*.


```r
x <- 5
x

y <- 42
y

y-x
z <- y-x
z

x <- 1:5
y <- 6:10
x
y
x+y
x*y
x^y
```

You can see what objects (variables) are stored by viewing the Environment tab in Rstudio. You can also use the `ls()` function. You can remove objects (variables) with the `rm()` function. You can do this one at a time or remove several objects at once.


```r
ls()
rm(x)
ls()
x # oops! you should get an error because x no longer exists!
```

```
## Error: object 'x' not found
```

```r
rm(y,z)
```

### Classes: everything is an object

Use the `class()` function to see what *kind* of object a something is. You can run `class()` on constants, built-in objects, or objects you create. Let's create a character object and then get it's class.


```r
class(42)
class(log)
name <- "Stephen"
name
class(name)
```

Certain *functions* operate only on certain *classes* of object. Here, `name` is a `character` class, assigned to `"Stephen"`. The built-in `toupper()` function will operate on character objects, but not others.


```r
toupper(name) # name is an object of class character. methods or functions are associated with certain classes.
toupper(log) # can't run a function that expects character on an object of class function
```

```
## Error: cannot coerce type 'special' to vector of type 'character'
```

We can combine values into a vector with the built-in `c()` function.


```r
# Get some help with ?c
x <- c(1,3,5)
x
class(x)
length(x)
```

Let's create and manipulate a character vector:


```r
y <- c("My", "name", "is", "Stephen")
y
class(y)
length(y)
y <- c(y, "Turner")
y
length(y)
```

Try running the built-in `sum()` function on a numeric vector:


```r
sum(x)
```

Combining characters with numerics results in coercing everything to be a character class.


```r
y
z <- c(x,y)
z
class(z)
```

Certain functions only operate on certain classes. You can't compute the `sum()` of a character vector!


```r
z
sum(z)
```

```
## Error: invalid 'type' (character) of argument
```

### Slicing/indexing vectors

Let's create a vector of 50 integers going from 101 to 150. We can access certain elements of that vector by putting the element's *index(es)* in square brackets. E.g., `x[1]` will return the first element in vector `x`. Calling `x[c(3,5)]` will access the third and fifth elements. Calling `x[1:10]` will return the first ten elements of `x`.

*Special note: R indexes vectors starting at 1. This is different from many other languages, including Perl and Python, which index starting from zero.*


```r
# Create the vector.
x <- 101:150

# Get the first element.
x[1]

# Get the 42nd element.
x[42]

# Get the 20th through the 25th elements.
x[20:25]

# If you try to access elements that don't exist, you'll return missing values.
# Missing values are represented as NA
x[45:55] #NA is missing value!
```

### Data Frames

Data frames are a standard way to store heterogeneous tabular data in R: tabular, meaning that individuals or observations are typically represented in rows, while variables or features are represented as columns; heterogeneous, meaning that columns/features/variables can be different classes (on variable, e.g. age, can be numeric, while another, e.g., cause of death, can be text).

Later on we'll go over how we load our own data, but for now, let's use a built-in data frame called `mtcars`. This data was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models). We can load this built-in data with `data(mtcars)`. By the way, running `data()` without any arguments will list all the available built-in datasets included with R.

Let's load the data first. Type the name of the object itself (`mtcars`) to view the entire data frame. *Note: doing this with large data frames can cause you trouble.*


```r
data(mtcars)
class(mtcars)
mtcars
```

There are several built-in functions that are useful for working with data frames.
* `head()` prints the first few lines of a large data frame.
* `length()` tells you the number of features (variables, columns) in a data frame.
* `dim()` returns a two-element vector containing the number of rows and the number of columns in a data frame.
* `str()` displays the structure of a data frame, printing out details and a preview of every column.
* `summary()` works differently depending on what kind of object you pass to it. Passing a data frame to the `summary()` function prints out some summary statistics about each column (min, max, median, mean, etc.)


```r
head(mtcars)
length(mtcars)
dim(mtcars)
dim(mtcars)[1] # number of rows (individual cars in the survey)
dim(mtcars)[2] # number of columns (number of variables measured)
str(mtcars)
```

We can access individual variables within a data frame using the `$` operator, e.g., `mydataframe$specificVariable`. Let's print out the number of cylinders for every car, and calculate the average miles per gallon for ever car in the dataset (using the built-in `mean()` function).


```r
# display the number of cylinders for each car.
mtcars$cyl
# first display MPG for all vehicles, then calculate the average.
mtcars$mpg
mean(mtcars$mpg)
```

We can also access certain rows or columns of a dataset by providing multiple indices using the syntax `mydataframe[rows, columns]`. Let's get the first 4 rows and the first two rows (MPG and # cylinders) from the dataset:


```r
head(mtcars)
mtcars[1:4, 1:2]
```

We can also use the `subset()` function to return a subset of the data frame that meets a specific condition. The first argument is the data frame you want to subset. The second argument is a condition you must satisfy. If you want to satisfy *all* of multiple conditions, you can use the "and" operator, `&`. The "or" operator `|` (the pipe character, usually shift-backslash) will return a subset that meet *any* of the conditions.

The commands below will:

0. Return only cars with 6 cylinder engines.
0. Return only cars with greater than 6 cylinders.
0. Return only the cars that get at least 20 miles per gallon or have a displacement volume of less than 100cc.
0. Return cars with 6 cylinder engines, but using the `select=` argument, only the MPG and displacement columns. Note the syntax there -- we're passing a vector of variables created with the `c()` function to the `select=` argument, which only returns certain columns.
0. Return cars that have greater than or equal to 6 cylinders *and* get at least 15 miles per gallon, but display only the MPG, cylinders, and qsec columns (qsec is the 1/4 mile time).

Try some subsetting on your own.


```r
subset(mtcars, cyl==6)
subset(mtcars, cyl>6)
subset(mtcars, mpg>=20 | disp<100)
subset(mtcars, cyl==6, select=c(mpg, disp))
subset(mtcars, cyl>=6 & mpg>=15, select=c(mpg, cyl, qsec))
```

The `with()` function is particularly helpful. Let's say you wanted to compute some (senseless) value by computing the MPG times the number of cylinders divided by the car's displacement. You could access the dataset's variables using the `$` notation, or you could use `with()` to temporarily *attach* the data frame, and call the variables directly. The first argument to `with()` is the name of the data frame, and the second argument is all the stuff you'd like to do with the particular features in that data frame.

Try typing the following commands:


```r
# Display the number of cylinders.
mtcars$cyl
with(mtcars, cyl)

# Compute the senseless value described above. Both return the same results.
mtcars$mpg * mtcars$cyl / mtcars$disp
with(mtcars, mpg*cyl/disp)
```

### Plotting

Plotting a single numeric variable goes down the rows and plots a value on the y-axis for each observation (index) in the data frame.


```r
plot(mtcars$mpg)
```

This isn't a very useful figure. More appropriate might be a histogram. We can try to let R decide how many breaks to insert in the histogram, or we can set that manually. We can also set the color of the bars.



```r
hist(mtcars$mpg)
hist(mtcars$mpg, breaks=10)
hist(mtcars$mpg, breaks=10, col="black")
```

We can create a scatterplot between two variables with `plot(varX, varY)`.


```r
# This would also work, but let's use with().
# plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))
```

There are hundreds of plotting parameters you can use to make your plot look exactly like you want. Let's use a solid-filled point instead of an open circle with the `pch=` argument, color the points red with the `col=` argument, give it a title by passing a character object to the `main=` argument, and change the x and y axis titles with the `xlab=` and `ylab=` arguments, respectively. Let's go through this one step at a time.


```r
with(mtcars, plot(disp, mpg, pch=16))
with(mtcars, plot(disp, mpg, pch=16,  col="red"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement",
                  ylab="Fuel Economy (MPG)", xlab="Displacement (cu. in.)"))
```

Notice how on that last line I broke the command up into two lines for better readability. I broke the command at the comma separating arguments, and indented the following line for readability.

On your own, try plotting horsepower vs displacement for vehicles with more than 4 cylinders. Give the graph a title and label the axes. Make the points solid (hint, `pch=16`) blue (hint, `col="blue"`) circles. Your plot should look something like this:



### Reading in / writing out data

First, lets create a small dataset consisting of only 8 cylinder cars.


```r
mtcars_8cyl <- subset(mtcars, cyl==8)
mtcars_8cyl
```

Next, check what your working directory is with `getwd()` with no arguments, and look up some help for `write.table()` and `write.csv()`.


```r
getwd()
help(write.table)
help(write.csv)
```

Using RStudio, go to the Session menu, and select the directory (folder) you want to work from under the "Set Working Directory" menu. You can also do this manually with the `setwd()` command.


```r
getwd()
setwd("~/Desktop/R")
```

Once you've set your working directory either using RStudio or on the command line, save the new reduced data frame to a comma-separated file called `cars8.csv` using the `write.csv()` function.


```r
write.csv(mtcars_8cyl, file="cars8.csv")
```

Data can be loaded using the Tools -- Import Dataset -- From text file menu in R studio. Or you can also load a dataset manually using `read.table()` or `read.csv()`. First, read the help on these functions:


```r
help(read.table)
help(read.csv)
```

Here let's remove the dataset, and re-import it into an object called cars8 from the file we just saved.


```r
rm(mtcars_8cyl)
mtcars_8cyl
cars8 <- read.table(file="cars8.csv", header=TRUE, sep=",", row.names=1)
cars8
rm(cars8)
cars8 <- read.csv(file="cars8.csv", header=TRUE, row.names=1)
cars8
```

## Using R and Bioconductor for gene expression analysis

In this section we'll analyze some publicly available gene expression data using R and bioinformatics-focused R packages in [Bioconductor](http://bioconductor.org/). But first, a bit about R *packages*.

### Packages

Most generic R packages are hosted on the Comprehensive R Archive Network (CRAN, <http://cran.us.r-project.org/>). To install one of these packages, you would use `install.packages("packagename")`. You only need to install a package once, then load it each time using `library(packagename)`. Let's install the **ggplot2** package, and load it.


```r
# Install only once.
install.packages("ggplot2")

# Load the package every time you want to use it.
library(ggplot2)
```

### Bioconductor

Bioconductor packages work a bit different, and are not hosted on CRAN. Go to <http://bioconductor.org/> to learn more about the Bioconductor project. To use any Bioconductor package, you'll need a few "core" Bioconductor packages. Run the following commands to (1) download the installer script, and (2) install some core Bioconductor packages. You'll need internet connectivity to do this, and it'll take a few minutes, but it only needs to be done once.


```r
# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function.
# Run it without any arguments to install the core packages or update any installed packages.
# This requires internet connectivity and will take some time!
biocLite()
```

To install specific packages, first download the installer script if you haven't done so, and use `biocLite("packagename")`. This only needs to be done once then you can load the package like any other package. Let's download the [DESeq2 package](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html):


```r
# Do only once
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

# Every time you need to use the DESeq2 package
library(DESeq2)
```

Bioconductor packages usually have great documentation in the form of *vignettes*. For a great example, take a look at the [DESeq2 vignette for analyzing count data](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).



### Analyzing publicly available RNA-seq data

Now, let's analyze some publicly available gene expression (RNA-seq) data. NCBI Gene Expression Omnibus (<http://www.ncbi.nlm.nih.gov/geo/>) is an international public repository that archives and freely distributes microarray, next-generation sequencing, and other forms of high-throughput functional genomics data submitted by the research community. Many publishers require gene expression data be submitted to GEO and made publicly available before publication. You can learn a lot more about GEO by reading their [overview](http://www.ncbi.nlm.nih.gov/geo/info/overview.html) and [FAQ](http://www.ncbi.nlm.nih.gov/geo/info/faq.html) pages. At the time of this writing, GEO hosts over 45,000 studies comprising over 1,000,000 samples on over 10,000 different technology platforms.

In this demonstration, we're going to be using data from GEO Series accession number GSE18508. You can enter this number in the search box on the GEO homepage, or use [this direct link](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). In this study, the authors performed RNA-sequencing of mRNA from *Drosophila melanogaster* S2-DRSC cells that have been RNAi depleted of mRNAs encoding RNA binding proteins and splicing factors. This was done as part of the modENCODE project, published in: Brooks AN et al. Conservation of an RNA regulatory map between Drosophila and mammals. *Genome Res* 2011 Feb;21(2):193-202.

You can go to the [GEO accession page](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508) and download all the raw data at the bottom under supplementary data. However, this dataset is nearly 100GB. What I have in this repository is a spreadsheet containing a matrix of gene counts, where genes are in rows and samples are in columns. That is, this data has already been aligned and counted, and the number in the cell is the number of RNA-seq reads that mapped to that gene for that sample. The value in the *i*-th row and the *j*-th column of the matrix tells how many reads have been mapped to gene *i* in sample *j*. To do these steps yourself, you would need to align reads to the genome (e.g., using [STAR](https://code.google.com/p/rna-star/)) and count reads mapping to features (e.g., using a GTF file from Ensembl and a tool like [featureCounts](http://bioinf.wehi.edu.au/featureCounts/))

Note: much of this was adapted from the [DESeq2 package vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf).

#### Load packages

First, we'll need to load the Bioconductor packages we'll be using:


```r
# Bioconductor packages.
# Use the installation instructions at http://www.bioconductor.org/install/
# if you haven't already installed these (install once, load every time)
library(Biobase)
library(DESeq2)
```

#### Load the data

Next, load two different datasets: one with the count data, and one with sample information. Take a look at the first few rows of the count dataset, as well as the entire metadata data frame.


```r
# Load the count data
pasillacounts <- read.csv("data/pasilla_counts.csv", header=TRUE, row.names=1)
head(pasillacounts)

# Load the sample metadata
pasillameta <- read.csv("data/pasilla_metadata.csv", header=TRUE, row.names=1)
pasillameta
```

The class used by the DESeq2 package to store the read counts is the **DESeqDataSet** class. This facilitates preparation steps andalso downstream exploration of results. A DESeqDataSet object must have an associated design formula. The design formula expresses the variables which will be used in modeling. The formula should be a tilde (~) followed by the variables with plus signs between them (it will be coerced into an formula if it is not already). An intercept is included, representing the base mean of counts. The design can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model. Do do this, we run the `DESeqDataSetFromMatrix` function, giving it data frames containing the count matrix as well as the column (metadata) information, and the design formula.


```r
dds <- DESeqDataSetFromMatrix(countData=pasillacounts, colData=pasillameta, design = ~condition)
dds
```

#### Differential expression analysis

The standard differential expression analysis steps are wrapped into a single function, `DESeq`. This convenience function normalizes the library, estimates dispersion, and performs a statistical test under the negative binomial model.


```r
# Normalize, estimate dispersion, fit model
dds <- DESeq(dds)
```

The results are accessed using the function `results`. The `mcols` function gives you more information about what the columns in the results tell you. You can reorder and write the results to a file if you wish. Finally, you can create an MA-plot, showing the log2 fold change over the normalized counts for each gene, colored red if statistically significant (padj<0.1).


```r
# Extract results
res <- results(dds)

# View, reorder, and write out results
head(res)
mcols(res)
res <- res[order(res$padj), ]
head(res)
write.csv(res, file="results/pasilla_results.csv")

# Create MA Plot
plotMA(dds, ylim=c(-2,2))
```

#### Data transformation and visualization

The differential expression analysis above operates on the raw (normalized) count data. But for visualizing or clustering data as you would with a microarray experiment, you ned to work with transformed versions of the data. First, use a *regularlized log* transofmration while re-estimating the dispersion ignoring any information you have about the samples (`blind=TRUE`). Perform a principal components analysis and hierarchical clustering.


```r
# Transform
rld <- rlogTransformation(dds, blind=TRUE)

# Principal components analysis
plotPCA(rld, intgroup=c("condition", "type"))

# Hierarchical clustering analysis
distrl <- dist(t(assay(rld)))
plot(hclust(distrl))
```

#### Multifactor design

See the [DESeq2 vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf) for more information on multifactor designs.

#### Record package and version info with `sessionInfo()`

The `sessionInfo()` prints version information about R and any attached packages. It's a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.


```r
sessionInfo()
```

```
## R version 3.0.2 (2013-09-25)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
##
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
##
## attached base packages:
## [1] parallel  methods   stats     graphics  grDevices utils     datasets
## [8] base
##
## other attached packages:
##  [1] DESeq2_1.2.10           RcppArmadillo_0.4.300.0
##  [3] Rcpp_0.11.1             GenomicRanges_1.14.4
##  [5] XVector_0.2.0           IRanges_1.20.7
##  [7] Biobase_2.22.0          BiocGenerics_0.8.0
##  [9] knitr_1.6               BiocInstaller_1.12.1
##
## loaded via a namespace (and not attached):
##  [1] annotate_1.40.1      AnnotationDbi_1.24.0 DBI_0.2-7
##  [4] evaluate_0.5.5       formatR_0.10         genefilter_1.44.0
##  [7] grid_3.0.2           lattice_0.20-29      locfit_1.5-9.1
## [10] RColorBrewer_1.0-5   RSQLite_0.11.4       splines_3.0.2
## [13] stats4_3.0.2         stringr_0.6.2        survival_2.37-7
## [16] tools_3.0.2          XML_3.95-0.2         xtable_1.7-3
```
