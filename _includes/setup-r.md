**Note:** R and RStudio are separate downloads and installations. **R** is the underlying statistical computing environment, but using R alone is no fun. **RStudio** is a graphical integrated development environment that makes using R much easier. You need R installed before you install RStudio.

0. **Download data.** First, [click this link to download and a zip file of this repository](https://github.com/bioconnector/workshops/archive/master.zip), then extract it somewhere on your computer where you can easily find it later (e.g. your Desktop). It has the data you'll need to go through these examples.
0. **Install R.** You'll need R version 3.0.0 or higher. Download and install R for [Windows](http://cran.r-project.org/bin/windows/base/) or [Mac OS X](http://cran.r-project.org/bin/macosx/) (download the latest R-3.x.x.pkg file for your appropriate version of OS X).
0. **Install RStudio.** Download and install the latest stable version of [RStudio Desktop](http://www.rstudio.com/products/rstudio/download/). Alternatively, download the [RStudio Desktop v0.99 preview release](http://www.rstudio.com/products/rstudio/download/preview/)  (the 0.99 preview version has many nice new features that are especially useful for this particular workshop).
0. **Install R packages.** Launch RStudio (RStudio, *not R itself*). Ensure that you have internet access, then enter the following commands into the **Console** panel (usually the lower-left panel, by default). Note that these commands are case-sensitive. At any point (especially if you've used R/Bioconductor in the past), R may ask you if you want to update any old packages by asking `Update all/some/none? [a/s/n]:`. If you see this, type `a` at the propt and hit `Enter` to update any old packages. If you're using a Windows machine you might get some errors about not having permission to modify the existing libraries -- don't worry about this message. You can avoid this error altogether by running RStudio as an administrator.

```r
# Install packages from CRAN
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("knitr")
install.packages("rmarkdown")

# Install Bioconductor packages
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq2")
```

You can check that you've installed everything correctly by closing and reopening RStudio and entering the following commands at the console window:

```r
library(dplyr)
library(ggplot2)
library(tidyr)
library(DESeq2)
```

These commands may produce some notes or other output, but as long as they work without an error message, you're good to go. If you get a message that says something like: `Error in library(packageName) : there is no package called 'packageName'`, then the required packages did not install correctly. Please do not hesitate to email me _prior to the course_ if you are still having difficulty.
