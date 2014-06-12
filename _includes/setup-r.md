You must bring a laptop with the necessary software installed to the course. Please install the software below *prior to the course* - we will not have time during the workshop to troubleshoot installation issues. Please email if you have any trouble.

0. **Download data.** [Click this link to download and a zip file of this repository](https://github.com/bioconnector/workshops/archive/master.zip), then extract it somewhere on your computer where you can easily find it later (e.g. your Desktop). It has the data you'll need to go through these examples.
0. **Install R.** You'll need R version 3.0.0 or higher. Download and install R for [Windows](http://cran.r-project.org/bin/windows/base/) or [Mac OS X](http://cran.r-project.org/bin/macosx/) (download the latest R-3.x.x.pkg file for your appropriate version of OS X).
0. **Install RStudio.** Download and install RStudio Desktop: <http://www.rstudio.com/ide/download/desktop>.
0. **Install R packages.** Launch RStudio (RStudio, *not R itself*). Ensure that you have internet access, then enter the following commands into the **Console** panel (usually the lower-left panel, by default). Note that these commands are case-sensitive. At any point (especially if you've used R/Bioconductor in the past), R may ask you if you want to update any old packages by asking `Update all/some/none? [a/s/n]:`. If you see this, type `a` at the propt and hit `Enter` to update any old packages.

```r
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biobase")
biocLite("DESeq2")
```
