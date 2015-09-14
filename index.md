---
layout: page
---

# Workshops

## Upcoming workshops

- [Workshop Series: Workshop Series: Data Analysis, Manipulation, and Visualization with R](ws-r-series/)
- [Reproducible Reporting: Generating Dynamic Documents with R+RStudio](ws-r-repres/)

## Past workshops

- [RNA-seq workshop](ws-rnaseq-1day/)
- [Introduction to R for Life Scientists](ws-r-intro/)
- [Advanced data manipulation with R and dplyr](ws-r-advanced-manipulation/)
- [Advanced Data Visualization with R and ggplot2](ws-r-advanced-visualization/)
- [Workshop Series: Workshop Series: Data Analysis, Manipulation, and Visualization with R](ws-r-series/)


# Updates

{% for post in site.posts %}
  * {{post.date | date: "%b %-d, %Y" }}: [{{post.title}}]({{ post.url | prepend: site.baseurl }})
{% endfor %}

[Subscribe to updates via RSS]({{ "/feed.xml" | prepend: site.baseurl }})
