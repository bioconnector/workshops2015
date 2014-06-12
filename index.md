---
layout: page
---

# Course Material

* [Introduction to R for Life Scientists](lessons/intro-r-lifesci/)
* RNA-seq workshop (*coming soon*)


# Updates

{% for post in site.posts %}
  * {{post.date | date: "%b %-d, %Y" }}: [{{post.title}}]({{ post.url | prepend: site.baseurl }})
{% endfor %}

[Subscribe to updates via RSS]({{ "/feed.xml" | prepend: site.baseurl }})
