---
layout: page
---

## Links to course material

Lorem ipsum.

## Updates

{% for post in site.posts %}
  * {{post.date | date: "%b %-d, %Y" }}: [{{post.title}}]({{ post.url | prepend: site.baseurl }})
{% endfor %}

<!--
  <h1>Posts</h1>

  <ul class="posts">
    {% for post in site.posts %}
      <li>
        <span class="post-date">{{ post.date | date: "%b %-d, %Y" }}</span>
        <a class="post-link" href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
      </li>
    {% endfor %}
  </ul>

  <p class="rss-subscribe">subscribe <a href="{{ "/feed.xml" | prepend: site.baseurl }}">via RSS</a></p> -->
