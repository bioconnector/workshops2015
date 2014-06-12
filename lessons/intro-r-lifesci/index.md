---
layout: page
---

# Introduction to R for Life Scientists

This workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce the R statistical computing environment, including basic instruction in data types, variables, array manipulation, functions, data frames, data import/export, visualization, and using packages. At the end of the workshop, participants will see a live demonstration of a real biomedical application - analysis of publicly available RNA-seq data. This will demo (1) how to acquire publicly accessible data from NCBI Gene Expression Omnibus, and (2) how to use Bioconductor packages to import, process, QC, analyze, and visualize the results of the analysis. By the end of the workshop, participants will be able to use R for basic data manipulation and visualization, and will know where to look for further help and instruction. Participants will also be exposed to analyzing RNA-seq data. An advanced follow-on course will go through a gene expression data analysis in detail.

**Date**: June 25, 2014  
**Time**: 2:00-5:00pm  
**Location (tentative)**: Health Sciences Library Carter Classroom (downstairs one floor and to the right).

**Pre-requisites**: You must bring a laptop to the course with the necessary software installed (see instructions below).

**Link to slides (*check back after course*)**

**[Link to course material](01-intro-r/)**

## Before coming

Before coming please install the software as instructed below, and take the pre-workshop survey (this should take ~60 seconds).

### Software setup

You must bring a laptop with the necessary software installed to the course. Please install the software below *prior to the course* - we will not have time during the workshop to troubleshoot installation issues. Please email if you have any trouble.

0. **Download data.** [Click this link to download and a zip file of this repository](https://github.com/stephenturner/teaching/archive/master.zip), then extract it somewhere on your computer where you can easily find it later (e.g. your Desktop). It has the data you'll need to go through these examples.
0. **Install R.** You'll need R version 3.0.0 or higher. Download and install R for [Windows](http://cran.r-project.org/bin/windows/base/) or [Mac OS X](http://cran.r-project.org/bin/macosx/) (download the R-3.x.x.pkg file for your appropriate version of OS X).
0. **Install RStudio.** Download and install RStudio Desktop: <http://www.rstudio.com/ide/download/desktop>.
0. **Install Bioconductor packages.** Launch RStudio (RStudio, *not R itself*). Ensure that you have internet access, then enter the following commands into the **Console** panel (usually the lower-left panel, by default). Note that these commands are case-sensitive.

```r
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biobase")
biocLite("DESeq2")
```

### Take the pre-workshop survey

[Click here to take the pre-workshop survey](https://docs.google.com/forms/d/1Ef4r-5yTOZO-rMGyjZ5M-wP3Q_j2WPtkp1HM_ksApnw/viewform). It should take about a minute to complete.
