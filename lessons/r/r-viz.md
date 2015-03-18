---
layout: page
---



# Advanced Data Visualization with R

This workshop will cover fundamental concepts for creating effective data visualization and will introduce tools and techniques for visualizing large, high-dimensional data using R. We will review fundamental concepts for visually displaying quantitative information, such as using series of small multiples, avoiding "chart-junk," and maximizing the data-ink ratio. After briefly covering data visualization using base R graphics, we will introduce the ggplot2 package for advanced high-dimensional visualization. We will cover the grammar of graphics (geoms, aesthetics, stats, and faceting), and using ggplot2 to create plots layer-by-layer. Upon completing this lesson, learners will be able to use ggplot2 to explore a high-dimensional dataset by faceting and scaling scatter plots in small multiples.

This is not a beginner course. This workshop requires basic familiarity with R: data types including vectors and data frames, importing/exporting data, and plotting. You can refresh your R knowledge with [DataCamp's Intro to R](https://www.datacamp.com/courses/introduction-to-r) or [TryR from CodeSchool](http://tryr.codeschool.com/). 

<!--
This course also requires some experience with data manipulation using dplyr functions. If you already know some R, you can learn the basics of dplyr by spending 10 minutes with [the dplyr vignette](http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html).
-->

_Attribution:_ This course material is modified from lesson material from Jenny Bryan's [Stat 545 course at UBC](http://stat545-ubc.github.io/), [Software Carpentry](http://software-carpentry.org/), and [Data Carpentry](http://datacarpentry.org/).

## Before coming

You'll need to bring a laptop to the course with the software installed as detailed below.

{% include setup-r.md %}

## Review

First, a little review about our data.

### The gapminder data

We use **data frames** to store heterogeneous tabular data in R: tabular, meaning that individuals or observations are typically represented in rows, while variables or features are represented as columns; heterogeneous, meaning that columns/features/variables can be different classes (on variable, e.g. age, can be numeric, while another, e.g., cause of death, can be text).

Before coming, you should have downloaded the gapminder data. If you [downloaded this entire lesson repository](https://github.com/bioconnector/workshops/archive/master.zip), once you extract it you'll find it in `workshops/lessons/r/data/gapminder.csv`. Alternatively you can download it directly from <http://bioconnector.org/data/>. This dataset is an excerpt from the [Gapminder](http://www.gapminder.org/) data, that's [already been cleaned up to a degree](https://github.com/jennybc/gapminder). This particular dataset has 1704 observations on six variables:

* `country` a categorical variable (aka "factor") 142 levels
* `continent`, a categorical variable with 5 levels
* `year`: going from 1952 to 2007 in increments of 5 years
* `pop`: population
* `gdpPercap`: GDP per capita
* `lifeExp`: life expectancy

Let's load the data first. There are two ways to do this. You can use RStudio's menus to select a file from your computer (tools, import dataset, from text file). But that's not reproducible. The best way to do this is to save the data and the script you're using to the same place, and read in the data in using `read.csv`. It's important to tell R that the file has a header, which tells R the names of the columns. We tell this to R with the `header=TRUE` argument. 

Once we've loaded it we can type the name of the object itself (`gm`) to view the entire data frame. Alternatively, we can use some of the functions below to look at bits and pieces of the data frame.


```r
# Read the data from file
gm <- read.csv("data/gapminder.csv", header=TRUE)

# Alternatively, read directly from the web:
# gm <- read.csv(url("http://bioconnector.org/data/gapminder.csv"), header=TRUE)
```

There are several built-in functions that are useful for working with data frames.

* Content:
    * `head()`: shows the first 6 rows
    * `tail()`: shows the last 6 rows
* Size:
    * `dim()`: returns a 2-element vector with the number of rows in the first element, and the number of columns as the second element (the dimensions of the object)
    * `nrow()`: returns the number of rows
    * `ncol()`: returns the number of columns
* Summary:
    * `colnames()` (or just `names()`): returns the column names
    * `str()`: structure of the object and information about the class, length and content of each column
    * `summary()`: works differently depending on what kind of object you pass to it. Passing a data frame to the `summary()` function prints out some summary statistics about each column (min, max, median, mean, etc.)
    

```r
head(gm)
tail(gm)
dim(gm)
nrow(gm)
ncol(gm)
colnames(gm)
str(gm)
summary(gm)
```

<!--

### dplyr review: verbs and chaining with **%>%**

First, let's turn gm into a tbl_df.


```r
library(dplyr)
gm <- tbl_df(gm)

gm
```

Recall the dplyr verbs and the chain or pipe operator `%>%`.

0. filter
0. select
0. mutate
0. arrange
0. summarize
0. group_by

#### filter()

The `filter()` function allows you to select out only particular rows that match a condition. 


```r
# Show only stats for the year 1982
filter(gm, year==1982)

# Alternatively,
gm %>% filter(year==1982)

# Show only stats for the US
gm %>% filter(country=="United States")

# Show only stats for American contries in 1997
gm %>% filter(continent=="Americas" & year==1997)
```

#### select()

The `filter()` function allows you to return only certain rows matching a condition. The `select()` function lets you subset the data and restrict to a number of columns. The first argument is the data, and subsequent arguments are the columns you want. Let's just get the year and the population variables.


```r
gm %>% select(year, pop)
```

#### mutate()

The `mutate()` function adds new columns to the data. Remember, the variable in our dataset is GDP per capita, which is the total GDP divided by the population size for that country, for that year. Let's mutate this dataset and add a column called gdpBil that is the raw GDP in billions:


```r
gm %>% mutate(gdpBil=pop*gdpPercap/1e9)
```


#### arrange()

The `arrange()` function does what it sounds like. It takes a data frame or tbl and arranges (or sorts) by column(s) of interest. The first argument is the data, and subsequent arguments are columns to sort on. Use the `desc()` function to arrange by descending.


```r
gm %>% arrange(lifeExp)
gm %>% arrange(year, desc(lifeExp))
```


#### group_by() %>% summarize()

Combining `group_by()` with `summarize()` can be very powerful. 


```r
gm %>% summarize(mean(lifeExp))
gm %>%
  group_by(continent, year) %>%
  summarize(mean(lifeExp))
```

--> 

## ggplot2 

**ggplot2** is a widely used R package that extends R's visualization capabilities. It takes the hassle out of things like creating legends, mapping other variables to scales like color, or faceting plots into small multiples. We'll learn about what all these things mean shortly. 

The **ggplot2** package provides an R implementation of Leland Wilkinson's *Grammar of Graphics* (1999). The *Grammar of Graphics* allows you to think beyond the garden variety plot types (e.g. scatterplot, barplot) and the consider the components that make up a plot or graphic, such as how data are represented on the plot (as lines, points, etc.), how variables are mapped to coordinates or plotting shape or color, what transformation or statistical summary is required, and so on. 

Specifically, **ggplot2** allows you to build a plot layer-by-layer by specifying:

- a **geom**, which specifies how the data are represented on the plot (points, lines, bars, etc.),
- **aesthetics** that map variables in the data to axes on the plot or to plotting size, shape, color, etc.,
- a **stat**, a statistical transformation or summary of the data applied prior to plotting,
- **facets**, which we've already seen above, that allow the data to be divided into chunks on the basis of other categorical or continuous variables and the same plot drawn for each chunk.

_First, a note about `qplot()`._ The `qplot()` function is a quick and dirty way of making ggplot2 plots. You might see it if you look for help with ggplot2, and it's even covered extensively in the ggplot2 book. And if you're used to making plots with built-in base graphics, the `qplot()` function will probably feel more familiar. But the sooner you abandon the `qplot()` syntax the sooner you'll start to really understand ggplot2's approach to building up plots layer by layer. So we're not going to use it at all in this class.

## Plotting bivariate data: continuous Y by continuous X

The `ggplot` function has two required arguments: the *data* used for creating the plot, and an *aesthetic* mapping to describe how variables in said data are mapped to things we can see on the plot.

First let's load the package:


```r
library(ggplot2)
```

Now, let's lay out the plot. If we want to plot a continuous Y variable by a continuous X variable we're probably most interested in a scatter plot. Here, we're telling ggplot that we want to use the `gm` dataset, and the aesthetic mapping will map `gdpPercap` onto the x-axis and `lifeExp` onto the y-axis.


```r
ggplot(gm, aes(x = gdpPercap, y = lifeExp))
```

```
## Error: No layers in plot
```

Look at that, we get an error, and it's pretty clear from the message what the problem is. We've laid out a two-dimensional plot specifying what goes on the x and y axes, but we haven't told it what kind of geometric object to plot. The obvious choice here is a point. Check out [docs.ggplot2.org](http://docs.ggplot2.org/) to see what kind of geoms are available.


```r
ggplot(gm, aes(x = gdpPercap, y = lifeExp)) + geom_point()
```

Here, we've built our plot in layers. First, we create a canvas for plotting layers to come using the `ggplot` function, specifying which **data** to use (here, the *gm* data frame), and an **aesthetic mapping** of `gdpPercap` to the x-axis and `lifeExp` to the y-axis. We next add a layer to the plot, specifying a **geom**, or a way of visually representing the aesthetic mapping. 

Now, the typical workflow for building up a ggplot2 plot is to first construct the figure and save that to a variable (for example, `p`), and as you're experimenting, you can continue to re-define the `p` object as you develop "keeper commands".

First, let's construct the graphic. Notice that we don't have to specify `x=` and `y=` if we specify the arguments in the correct order (x is first, y is second).


```r
p <- ggplot(gm, aes(gdpPercap, lifeExp))
```

Now, if we tried to display p here alone we'd get another error because we don't have any layers in the plot. Let's experiment with adding points and a different scale to the x-axis.


```r
# Experiment with adding poings
p + geom_point()

# Experiment with a different scale
p + geom_point() + scale_x_log10()
```

I like the look of using a log scale for the x-axis. Let's make that stick.


```r
p <- p + scale_x_log10()
```

Then re-plot again with a layer of points:


```r
p + geom_point()
```

Now notice what I've saved to `p` at this point: only the basic plot layout and the log10 mapping on the x-axis. I didn't save any layers yet because I want to fiddle around with the points for a bit first.

Above we implied the aesthetic mappings for the x- and y- axis should be `gdpPercap` and `lifeExp`, but we can also add aesthetic mappings to the geoms themselves. For instance, what if we wanted to color the points by the value of another variable in the dataset, say, continent?


```r
p + geom_point(aes(color=continent))
```

Notice the difference here. If I wanted the colors to be some static value, I wouldn't wrap that in a call to `aes()`. I would just specify it outright. Same thing with other features of the points. For example, lets make all the points huge (`size=8`) blue (`color="blue"`) semitransparent (`alpha=(1/4)`) triangles (`pch=17`):


```r
p + geom_point(color="blue", pch=17, size=8, alpha=1/4)
```

Now, this time, let's map the aesthetics of the point character to certain features of the data. For instance, let's give the points different colors and character shapes according to the continent, and map the size of the point onto the life Expectancy:


```r
p + geom_point(aes(col=continent, pch=continent, size=lifeExp))
```

Now, this isn't a great plot because there are several aesthetic mappings that are redundant. Life expectancy is mapped to both the y-axis and the size of the points -- the size mapping is superfluous. Similarly, continent is mapped to both the color and the point character (the point character is superfluous). Let's get rid of that, but let's make the points a little bigger outsize of an aesthetic mapping.


```r
p + geom_point(aes(col=continent), size=4)
```

---

**EXERCISE 1**

Re-create this same plot from scratch without saving anything to a variable. That is, start from the `ggplot` call. 

* Start with the `ggplot()` function.
* Use the gm data.
* Map `gdpPercap` to the x-axis and `lifeExp` to the y-axis.
* Add points to the plot
  * Make the points size 4
  * Map continent onto the aesthetics of the point
* Use a log<sub>10</sub> scale for the x-axis.



---

### Adding layers to the plot

Let's add a fitted curve to the points. 


```r
p <- ggplot(gm, aes(gdpPercap, lifeExp)) + scale_x_log10()
p + geom_point() + geom_smooth()
```

By default `geom_smooth()` will try to lowess for data with n<1000 or generalized additive models for data with n>1000. We can change that behavior by tweaking the parameters to use a thick red line, use a linear model instead of a GAM, and to turn off the standard error stripes.


```r
p + geom_point() + geom_smooth(lwd=2, se=FALSE, method="lm", col="red")
```

But let's add back in our aesthetic mapping to the continents. Notice what happens here. We're mapping continent as an aesthetic mapping _to the color of the points only_ -- so `geom_smooth()` still works only on the entire data. 


```r
p + geom_point(aes(color = continent)) + geom_smooth()
```

But notice what happens here: we make the call to `aes()` outside of the `geom_point()` call, and the continent variable gets mapped as an aesthetic to any further geoms. So here, we get separate smoothing lines for each continent. Let's do it again but remove the standard error stripes and make the lines a bit thicker.


```r
p + aes(color = continent) + geom_point() + geom_smooth()
p + aes(color = continent) + geom_point() + geom_smooth(se=F, lwd=2)
```

### Faceting

Facets display subsets of the data in different panels. There are a couple ways to do this, but `facet_wrap()` tries to sensibly wrap a series of facets into a 2-dimensional grid of small multiples. Just give it a formula specifying which variables to facet by. We can continue adding more layers, such as smoothing. If you have a look at the help for `?facet_wrap()` you'll see that we can control how the wrapping is laid out.


```r
p + geom_point() + facet_wrap(~continent)
p + geom_point() + geom_smooth() + facet_wrap(~continent)
p + geom_point() + geom_smooth() + facet_wrap(~continent, ncol=1)
```

### Saving plots

There are a few ways to save ggplots. The quickest way, that works in an interactive session, is to use the `ggsave()` function. You give it a file name and by default it saves the last plot that was printed to the screen. 


```r
p + geom_point()
ggsave(file="myplot.png")
```

But if you're running this through a script, the best way to do it is to pass `ggsave()` the object containing the plot that is meant to be saved. We can also adjust things like the width, height, and resolution. `ggsave()` also recognizes the name of the file extension and saves the appropriate kind of file. Let's save a PDF.


```r
pfinal <- p + geom_point() + geom_smooth() + facet_wrap(~continent, ncol=1)
ggsave(pfinal, file="myplot.pdf", width=5, height=15)
```

---

**EXERCISE 2**

0. Make a scatter plot of `lifeExp` on the y-axis against `year` on the x.
0. Make a series of small multiples faceting on continent.
0. Add a fitted curve, smooth or lm, with and without facets.
0. **Bonus**: using `geom_line()` and and aesthetic mapping `country` to `group=`, make a "spaghetti plot", showing _semitransparent_ lines connected for each country, faceted by continent. Add a smoothed loess curve with a thick (`lwd=3`) line with no standard error stripe. Reduce the opacity (`alpha=`) of the individual black lines.



---

## Plotting bivariate data: continuous Y by categorical X

With the last example we examined the relationship between a continuous Y variable against a continuous X variable. A scatter plot was the obvious kind of data visualization. But what if we wanted to visualize a continuous Y variable against a categorical X variable? We sort of saw what that looked like in the last exercise. `year` is a continuous variable, but in this dataset, it's broken up into 5-year segments, so you could almost think of each year as a categorical variable. But a better example would be life expectancy against continent or country. 

First, let's set up the basic plot:


```r
p <- ggplot(gm, aes(continent, lifeExp)) 
```

Then add points:


```r
p + geom_point()
```

That's not terribly useful. There's a big overplotting problem. We can try to solve with transparency:


```r
p + geom_point(alpha=1/4)
```

But that really only gets us so far. What if we spread things out by adding a little bit of horizontal noise (aka "jitter") to the data.


```r
p + geom_jitter()
```

Note that the little bit of horizontal noise that's added to the jitter is random. If you run that command over and over again, each time it will look slightly different. The idea is to visualize the density at each vertical position, and spreading out the points horizontally allows you to do that. If there were still lots of over-plotting you might think about adding some transparency by setting the `alpha=` value for the jitter.


```r
p + geom_jitter(alpha=1/2)
```

Probably a more common visualization is to show a box plot:


```r
p + geom_boxplot()
```

But why not show the summary and the raw data?


```r
p + geom_jitter() + geom_boxplot()
```

Notice how in that example we first added the jitter layer then added the boxplot layer. But the boxplot is now superimposed over the jitter layer. Let's make the jitter layer go on top. Also, go back to just the boxplots. Notice that the outliers are represented as points. But there's no distinction between the outlier point from the boxplot geom and all the other points from the jitter geom. Let's change that. Notice the British spelling.


```r
p + geom_boxplot(outlier.colour = "red") + geom_jitter(alpha=1/2)
```

There's another geom that's useful here, called a voilin plot.


```r
p + geom_violin()

p + geom_violin() + geom_jitter(alpha=1/2)
```

Let's go back to our boxplot for a moment.


```r
p + geom_boxplot()
```

This plot would be a lot more effective if the continents were shown in some sort of order other than alphabetical. To do that, we'll have to go back to our basic build of the plot again and use the `reorder` function in our original aesthetic mapping. Here, reorder is taking the first variable, which is some categorical variable, and ordering it by the level of the mean of the second variable, which is a continuous variable. It looks like this


```r
p <- ggplot(gm, aes(x=reorder(continent, lifeExp), y=lifeExp))
```


```r
p + geom_boxplot()
```

---

**EXERCISE 3**

0. Make a jittered strip plot of GDP per capita against continent.
0. Make a box plot of GDP per capita against continent.
0. Using a log<sub>10</sub> y-axis scale, overlay semitransparent jittered points on top of box plots, where outlying points are colored. 
0. **BONUS**: Try to reorder the continents on the x-axis by GDP per capita. Why isn't this working as expected? See `?reorder` for clues.



---

## Plotting univariate continuous data

What if we just wanted to visualize distribution of a single continuous variable? A histogram is the usual go-to visualization. Here we only have one aesthetic mapping instead of two.


```r
p <- ggplot(gm, aes(lifeExp))
p + geom_histogram()
```

When we do this ggplot lets us know that we're automatically selecting the width of the bins, and we might want to think about this a little further.


```r
p + geom_histogram(binwidth=5)
p + geom_histogram(binwidth=1)
p + geom_histogram(binwidth=.25)
```

Alternative we could plot a smoothed density curve instead of a histogram:


```r
p + geom_density()
```

Back to histograms. What if we wanted to color this by continent?


```r
p + geom_histogram(aes(color=continent))
```

That's not what we had in mind. That's just the outline of the bars. We want to change the fill color of the bars.


```r
p + geom_histogram(aes(fill=continent))
```

Well, that's not exactly what we want either. If you look at the help for `?geom_histogram` you'll see that by default it stacks overlapping points. This isn't really an effective visualization. Let's change the position argument.


```r
p + geom_histogram(aes(fill=continent), position="identity")
```

But the problem there is that the histograms are blocking each other. What if we tried transparency?


```r
p + geom_histogram(aes(fill=continent), position="identity", alpha=1/3)
```

That's somewhat helpful, and might work for two distributions, but it gets cumbersome with 5. Let's go back and try this with density plots, first changing the color of the line:


```r
p + geom_density(aes(color=continent))
```

Then by changing the color of the fill and setting the transparency to 25%:


```r
p + geom_density(aes(fill=continent), alpha=1/4)
```


---

**EXERCISE 4**

0. Plot a histogram of GDP Per Capita.
0. Do the same but use a log<sub>10</sub> x-axis.
0. Still on the log<sub>10</sub> x-axis scale, try a density plot mapping continent to the fill of each density distribution, and reduce the opacity.
0. Still on the log<sub>10</sub> x-axis scale, make a histogram faceted by continent _and_ filled by continent. Facet with a single column (see `?facet_wrap` for help). Save this to a 6x10 PDF file.




---

## Themes

Let's make a plot we made earlier (life expectancy versus the log of GDP per capita with points colored by continent with lowess smooth curves overlaid without the standard error ribbon):


```r
p <- ggplot(gm, aes(gdpPercap, lifeExp)) 
p <- p + scale_x_log10()
p <- p + aes(col=continent) + geom_point() + geom_smooth(lwd=2, se=FALSE)
p
```

Give the plot a title and axis labels:


```r
p <- p + ggtitle("Life expectancy vs GDP by Continent")
p <- p + xlab("GDP Per Capita (USD)") + ylab("Life Expectancy (years)")
p
```

They "gray" theme is the usual background.


```r
p + theme_gray()
```

We could also get a black and white background:


```r
p + theme_bw()
```

Or go a step further and remove the gridlines:


```r
p + theme_classic()
```

Finally, there's another package that gives us lots of different themes. This package isn't on CRAN, so you'll have to use devtools to install it directly from the source code on github.


```r
install.packages("devtools")
devtools::install_github("jrnold/ggthemes")
```


```r
library(ggthemes)
p + theme_excel()
p + theme_excel() + scale_colour_excel()
p + theme_gdocs() + scale_colour_gdocs()
p + theme_stata() + scale_colour_stata()
p + theme_wsj() + scale_colour_wsj()
p + theme_economist() 
p + theme_fivethirtyeight()
p + theme_tufte()
```



<!-- This begins a comment

---

**EXERCISE 5**



---

comment ends here -->


## Further ggplot2 resources

* <http://docs.ggplot2.org/>: The official ggplot2 documentation.
* <http://amzn.to/1akjqsR>: Edition 1 of the ggplot2 book, by the developer, Hadley Wickham.
* <https://github.com/hadley/ggplot2-book>: New version of the ggplot2 book, freely available on GitHub.
* <https://groups.google.com/d/forum/ggplot2>: The ggplot2 Google Group (mailing list, discussion forum).
* <http://learnr.wordpress.com/>: A blog with a good number of posts describing how to reproduce various kind of plots using ggplot2.
* <http://stackoverflow.com/questions/tagged/ggplot2>: Thousands of questions and answers tagged with "ggplot2" on Stack Overflow, a programming Q&A site.
* <http://stat545-ubc.github.io/>: Jenny Bryan's stat 545 class at UBC -- a great resource for learning about data manipulation and visualization in R. Much of this course material was modified from the stat 545 lesson material.
* <http://shinyapps.stat.ubc.ca/r-graph-catalog/>: A catalog of graphs made with ggplot2, complete with accompanying R code.
