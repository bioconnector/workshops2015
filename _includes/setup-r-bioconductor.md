Additionally, you'll need to install a few [Bioconductor](http://bioconductor.org/) packages. These packages are installed differently than "regular" R packages from CRAN. Copy and paste these lines of code into your R console.

```r
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq2")
```

You can check that you've installed everything correctly by closing and reopening RStudio and entering the following commands at the console window:

```r
library(DESeq2)
```

If you get a message that says something like: `Error in library(packageName) : there is no package called 'packageName'`, then the required packages did not install correctly. Please do not hesitate to email me _prior to the course_ if you are still having difficulty.
