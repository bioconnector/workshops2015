---
layout: page
---

## Link to Course Material

[Click here to access all course material](lessons).


## Updates

{% for post in site.posts %}
  * {{post.date | date: "%b %-d, %Y" }}: [{{post.title}}]({{ post.url | prepend: site.baseurl }})
{% endfor %}

[Subscribe to updates via RSS]({{ "/feed.xml" | prepend: site.baseurl }})
