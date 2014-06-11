---
layout: page
root: ../..
---

# Introduction to R for Life Scientists

This workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce the R statistical computing environment, including basic instruction in data types, variables, array manipulation, functions, data frames, data import/export, visualization, and using packages. At the end of the workshop, participants will see a live demonstration of a real biomedical application - analysis of publicly available RNA-seq data. This will demo (1) how to acquire publicly accessible data from NCBI Gene Expression Omnibus, and (2) how to use Bioconductor packages to import, process, QC, analyze, and visualize the results of the analysis. By the end of the workshop, participants will be able to use R for basic data manipulation and visualization, and will know where to look for further help and instruction. Participants will also be exposed to analyzing RNA-seq data. An advanced follow-on course will go through a gene expression data analysis in detail.

**Date**: June 25, 2014  
**Time**: 2:00-5:00pm  
**Location**: Health Sciences Library Carter Classroom (downstairs one floor and to the right).

**Registration**: [See below](#how-to-register).

**Pre-requisites**: You must bring a laptop to the course with the necessary software installed (see instructions below). You must install this software before you can register.

**Link to slides (*check back after course*)**

**[Link to course material](01-intro-r/)**


## How to register

0. Install the software and run the code as detailed below. This will reveal a link to take the pre-workshop survey.
0. Go to the link to take the pre-workshop survey (should take you ~60 seconds). Completing the survey will reveal the actual registration link.
0. Go to this link and register.

Please don't hesitate to [email me](http://stephenturner.us/email) if you have trouble getting set up. This process will ensure that everyone is ready to go so we don't have to spend the first 20 minutes getting people set up.

### Software setup

0. **Download data.** [Click this link to download and a zip file of this repository](https://github.com/stephenturner/teaching/archive/master.zip), then extract it somewhere on your computer where you can easily find it later (e.g. your Desktop). It has the data you'll need to go through these examples.
0. **Install R.** You'll need R version 3.0.0 or higher. Download and install R for [Windows](http://cran.r-project.org/bin/windows/base/) or [Mac OS X](http://cran.r-project.org/bin/macosx/) (download the R-3.x.x.pkg file for your appropriate version of OS X).
0. **Install RStudio.** Download and install RStudio Desktop: <http://www.rstudio.com/ide/download/desktop>.
0. **Install packages.** Launch RStudio (RStudio, *not R itself*). Ensure that you have internet access, then enter the following commands into the **Console** panel (usually the lower-left panel, by default). Note that these commands are case-sensitive.
  * The first 4 commands will download and install necessary add-on packages we will use in class. If anything asks you to update, type "**a**" and hit enter at the prompt to update **a**ll packages.
  * The second 4 commands will install a package I wrote, and runs a function that will display the link to take the survey.

```coffee
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biobase")
biocLite("DESeq2")

install.packages("devtools")
library(devtools)
install_github("stephenturner/Tmisc")
library(Tmisc)
registration()
```

### Take the pre-workshop survey

Successfully running the commands above will reveal a link to take a pre-workshop survey. This survey should take about one minute to complete. After you submit the survey you'll get a link to register.

Again, please don't hesitate to [email me](http://stephenturner.us/email) if you have trouble getting set up.
