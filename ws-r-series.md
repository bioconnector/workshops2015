---
layout: page
---

# Workshop Series: Data Analysis, Manipulation, and Visualization with R

This is a three-part series. The first session will introduce participants to the R environment and a dataset to be used for the remainder of the series. The two latter sessions will introduce participants to advanced data manipulation and visualization, using the same dataset introduced in the first session. This three-part workshop series will use a consistent example dataset for answering questions and completing exercises across all three sessions, culminating with a "capstone" analysis that integrates all covered material.


**Pre-requisites**:

- Bring a laptop to the course with the software installed [as detailed below](#setup).
- Print these cheat sheets and bring them with you to the workshop:
    - [R Reference Card](http://cran.r-project.org/doc/contrib/Short-refcard.pdf)
    - [Advanced Data Wrangling Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
    - [Advanced Data Visualization Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2015/04/ggplot2-cheatsheet.pdf)

**Registration**: [See below](#registration) for registration instructions.

**Instructor / Technical contact**: Stephen Turner  (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">s...</a>@virginia.edu)  
**Logistics / registration contact**: Bart Ragon (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75Vsnuy3VwvY13wVeE0K2DFU5Cf-2n-YnO3260iwqa1RA\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">b...</a>@virginia.edu)

## Syllabus

This workshop is a three-part series. Each subsequent workshop builds on the material covered in earlier sessions.

**Part I:** Sepbember 21, 1:00pm - 5:00pm  
**Part II:** Sepbember 22, 1:00pm - 5:00pm  
**Part II:** Sepbember 23, 1:00pm - 5:00pm

**Location**: Carter classroom, first floor Health Sciences Library (down the stairs and to the right).

***Instruction will start promptly at 1:00pm on the first day***. If you have any trouble with setup, please contact Stephen Turner prior to the course.

### Part I: Introduction to R

This beginner-level workshop is directed toward life scientists with little to no experience with statistical computing or bioinformatics. This interactive workshop will introduce the R statistical computing environment, The first part of this workshop will demonstrate very basic functionality in R, including functions, functions, vectors, creating variables, getting help, filtering, data frames, plotting, and reading/writing files.

### Part II: Advanced Data Manipulation with R

Data analysis involves a large amount of [janitor work](http://www.nytimes.com/2014/08/18/technology/for-big-data-scientists-hurdle-to-insights-is-janitor-work.html) - munging and cleaning data to facilitate downstream data analysis. This workshop is designed for those with a basic familiarity with R who want to learn tools and techniques for advanced data manipulation. It will cover data cleaning and "tidy data," and will introduce participants to R packages that enable data manipulation, analysis, and visualization using split-apply-combine strategies. Upon completing this lesson, participants will be able to use the dplyr package in R to effectively manipulate and conditionally compute summary statistics over subsets of a "big" dataset containing many observations.

### Part III: Advanced Data Visualization with R and ggplot2

This workshop will cover fundamental concepts for creating effective data visualization and will introduce tools and techniques for visualizing large, high-dimensional data using R. We will review fundamental concepts for visually displaying quantitative information, such as using series of small multiples, avoiding "chart-junk," and maximizing the data-ink ratio. After briefly covering data visualization using base R graphics, we will introduce the ggplot2 package for advanced high-dimensional visualization. We will cover the grammar of graphics (geoms, aesthetics, stats, and faceting), and using ggplot2 to create plots layer-by-layer. Upon completing this lesson, learners will be able to use ggplot2 to explore a high-dimensional dataset by faceting and scaling scatter plots in small multiples.

## Course Material

_Check back after course_.

1. [Introduction to R](../lessons/r/r-intro/)

<!--
1. [Introduction to R](../lessons/r/r-intro/)
1. [Advanced Data Manipulation](../lessons/r/r-manipulation/)
1. [Advanced Data Visualization](../lessons/r/r-viz/)
-->

<a name="setup"></a>

## Setup

You must bring a laptop with the necessary software installed to the course. Please install the software below *prior to the course* - we will not have time during the workshop to troubleshoot installation issues. Please email me (<a href="http://www.google.com/recaptcha/mailhide/d?k=01uXi4zl-bIdygzSeXF4649A==&amp;c=_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE=" onclick="window.open('http://www.google.com/recaptcha/mailhide/d?k\07501uXi4zl-bIdygzSeXF4649A\75\75\46c\75_81hv-sTQvJ9rjELjZNDJeAXTvLvkpfD9KEuItpEHTE\075', '', 'toolbar=0,scrollbars=0,location=0,statusbar=0,menubar=0,resizable=0,width=500,height=300'); return false;" title="Reveal this e-mail address">sd...</a>@virginia.edu) if you have any trouble.

{% include setup-r.md %}

<a name="registration"></a>

## Registration

**[Click here to register](https://www.bioconnector.virginia.edu/workshops/workshop-series-data-analysis-manipulation-and-visualization-r).**

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
