---
layout: page
---

# RNA-seq data analysis bootcamp

This workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce both the Linux/UNIX operating system and the R statistical computing environment, with a focus on a biological application - analyzing RNA-seq data for differentially expressed genes. The morning session will introduce basic operation in a UNIX environment, and will cover the first steps in an RNA-seq analysis including QC, alignment, and quantitation. The afternoon will introduce the R statistical computing environment, and will cover differential gene expression analysis using Bioconductor. By the end of the workshop, participants will:

0. Be familiar with the UNIX shell, including nagivating the filesystem, creating/examining/removing files, getting help, and batch operations.
0. Know how to align and quantitate gene expression with RNA-seq data
0. Become familiar with the R statistical computing environment, including data types, variables, array manipulation, functions, data frames, data import/export, visualization, and using packages.
0. Know what packages to use and what steps to take to analyze RNA-seq data for differentially expressed genes.

Participants will also be exposed to operating in a virtual environment and/or provisioning their own cloud computing resources. This course is sponsored by the Claude Moore Health Sciences Library, and borrows some materials from the Software Carpentry and Data Carpentry projects.

**Date**: Monday, November 10, 2014  
**Time**: 8:00 am (sharp!) - 5:00 pm  
**Location**: Carter classroom, first floor Health Sciences Library

**Pre-requisites**: [See the setup requirements below](#setup). Set aside an hour to create the necessary accounts and install the software *prior to* the workshop. *We will not have time to do this during the workshop*.

**Instructor / Technical contact**: Stephen Turner  (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">s...</a>@virginia.edu)  
**Logistics / registration contact**: Bart Ragon (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">b...</a>@virginia.edu)

## Agenda

***The workshop will start promptly at 8:00am***. If you have any trouble with setup, please contact Stephen Turner prior to the course. Dr. Turner will also be available 30 minutes prior to the course for hands-on troubleshooting, but please try to solve any setup problems prior to this time if possible.

* 0730-0800: *(Optional)* Help with setup
* 0800-0945: Introduction to Linux
* 1000-1200: QC, alignment and expression quantitation
* 1200-1300: Lunch (provided)
* 1300-1445: Introduction to R
* 1500-1700: QC, differential expression, and visualization with R/Bioconductor

## Course Material

* [Introduction to Unix](../shell/01-intro-unix-shell/)
* [NGS data analysis: QC, Alignment, Quantitation](01-alignment-counting/)
* [R Slides](https://speakerdeck.com/stephenturner/introduction-to-r-for-life-scientists)
* [Introduction to R](../intro-r-lifesci/01-intro-r/)
* [RNA-seq data analysis with DESeq2](03-differential-expression/)

## Registration

Registration is $20 (includes lunch). Click the link below and scroll down to find the "Register" button, or use `Ctrl-F` in your browser and search for `RNA-seq Data Analysis Bootcamp`.

<http://www.medicine.virginia.edu/education/more/cme/all-cme-prog-folder/FLP-page>

<!--
    This block includes the Eventbrite registration widget if 'eventbrite' has been set in the header.

    Maybe you need to change height value:

    - for one room use 206px,
    - for one waitlist room use 152px,
    - for two room use 254px,
    - for one waitlist room and one room use 253px,
    - for two waitlist room use 197px.
-->

{% if page.eventbrite %}
<iframe src="//www.eventbrite.com/tickets-external?eid={{page.eventbrite}}&ref=etckt" frameborder="0" width="100%" height="253px" scrolling="auto"></iframe>
{% endif %}

<a name="setup"></a>

## Setup

Please bring a laptop with the software below installed (everything is free). You'll also need to create an Amazon Web Services account. I can't understate how important it is to do this *prior to the course* - we will not have time during the workshop to troubleshoot installation issues. Please email me (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">sd...</a>@virginia.edu) if you have any trouble.

Setup checklist:

* Register & activate an AWS account
* Get a free AWS voucher from Dr. Turner
* Download PuTTY (Windows users only)
* Install Cyberduck
* Download and extract course repository zip file
* Install R
* Install RStudio
* Install DESeq R package


### Software setup, part I: AWS and a terminal

{% include setup-shell.md %}

### Software setup, part II: R and RStudio

{% include setup-r.md %}
