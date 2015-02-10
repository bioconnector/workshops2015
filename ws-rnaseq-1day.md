---
layout: page
eventbrite: 15712713148
---

# RNA-seq data analysis bootcamp

This workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce both the Linux/UNIX operating system and the R statistical computing environment, with a focus on a biological application - analyzing RNA-seq data for differentially expressed genes. The morning session will introduce basic operation in a UNIX environment, and will cover the first steps in an RNA-seq analysis including QC, alignment, and quantitation. The afternoon will introduce the R statistical computing environment, and will cover differential gene expression analysis using Bioconductor. By the end of the workshop, participants will:

0. Be familiar with the UNIX shell, including nagivating the filesystem, creating/examining/removing files, getting help, and batch operations.
0. Know how to align and quantitate gene expression with RNA-seq data
0. Become familiar with the R statistical computing environment, including data types, variables, array manipulation, functions, data frames, data import/export, visualization, and using packages.
0. Know what packages to use and what steps to take to analyze RNA-seq data for differentially expressed genes.

Participants will also be exposed to operating in a virtual environment and/or provisioning their own cloud computing resources. This course is sponsored by the Claude Moore Health Sciences Library, and borrows some materials from the Software Carpentry and Data Carpentry projects.

**Pre-requisites**: [See the setup requirements below](#setup). Set aside an hour to create the necessary accounts and install the software *prior to* the workshop. *We will not have time to do this during the workshop*.

**Registration**: Registration opens Monday February 23, 2015 at 9:00am. [See below](#registration) for registration instructions.

**Instructor / Technical contact**: Stephen Turner  (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">s...</a>@virginia.edu)  
**Logistics / registration contact**: Bart Ragon (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">b...</a>@virginia.edu)

## Agenda

The boot camp is a two-part series.

**Part I:** Monday, March 23 2015, 8:30am - 12:30pm  
**Part II:** Thursday, March 26 2015, 1:00pm - 5:00pm  
**Location**: Carter classroom, first floor Health Sciences Library

***Instruction will start promptly at 8:30am on the first day***. If you have any trouble with setup, please contact Stephen Turner prior to the course. Dr. Turner will also be available at 8:00am that morning, 30 minutes prior to the course for hands-on troubleshooting, but please try to solve any setup problems prior to this time if possible.

**Part I:**

* 0800-0830: *(Optional)* Help with setup
* 0830-1030: Introduction to Linux
* 1045-1230: QC, alignment and expression quantitation

**Part II:**

* 1300-1445: Introduction to R
* 1500-1700: QC, differential expression, and visualization with R/Bioconductor

## Course Material

_Check back after course_.

<!--
* [Introduction to Unix](../lessons/shell/shell-intro/)
* [NGS data analysis: QC, Alignment, Quantitation](../lessons/rnaseq/rnaseq-align-count/)
* [R Slides](https://speakerdeck.com/stephenturner/introduction-to-r-for-life-scientists)
* [Introduction to R](../lessons/r/r-intro/)
* [RNA-seq data analysis with DESeq2](../lessons/rnaseq/rnaseq-diff-expr/)
 -->

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
* Install DESeq2 R package

### Software setup, part I: AWS and a terminal

{% include setup-shell.md %}

### Software setup, part II: R and RStudio

{% include setup-r.md %}

<a name="registration"></a>

## Registration

[Register here](http://www.eventbrite.com/e/rna-seq-boot-camp-tickets-15712713148) or use the form below. The registration fee is $10, and that covers both parts. Registration is non-refundable, and part I and II cannot be split and attended individually.

Registration opens Monday February 23, 2015 at 9:00am.

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
<iframe src="//www.eventbrite.com/tickets-external?eid={{page.eventbrite}}&ref=etckt" frameborder="0" width="100%" height="250px" scrolling="auto"></iframe>
{% endif %}
