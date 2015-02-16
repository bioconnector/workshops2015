---
layout: page
---



# Advanced Data Manipulation with R

Data analysis involves a large amount of [janitor work](http://www.nytimes.com/2014/08/18/technology/for-big-data-scientists-hurdle-to-insights-is-janitor-work.html) - munging and cleaning data to facilitate downstream data analysis. This workshop assumes you have a basic familiarity with R and want to learn tools and techniques for advanced data manipulation. It will cover data cleaning and "tidy data," and will introduce R packages that enable data manipulation, analysis, and visualization using split-apply-combine strategies. Upon completing this lesson, you will be able to use the dplyr package in R to effectively manipulate and conditionally compute summary statistics over subsets of a "big" dataset containing many observations.

This is not a beginner course. This workshop requires basic familiarity with R: data types including vectors and data frames, importing/exporting data, and plotting. You can refresh your R knowledge with [DataCamp's Intro to R](https://www.datacamp.com/courses/introduction-to-r) or [TryR from CodeSchool](http://tryr.codeschool.com/).

## Before coming

You'll need to bring a laptop to the course with the software installed as detailed below.

{% include setup-r.md %}

## Review

We use **data frames** to store heterogeneous tabular data in R: tabular, meaning that individuals or observations are typically represented in rows, while variables or features are represented as columns; heterogeneous, meaning that columns/features/variables can be different classes (on variable, e.g. age, can be numeric, while another, e.g., cause of death, can be text).

Before coming, you should have downloaded the gapminder data. If you [downloaded this entire lesson repository](https://github.com/bioconnector/workshops/archive/master.zip), once you extract it you'll find it in `workshops/lessons/r/data/gapminder.csv`. Alternatively you can download it directly from <http://bioconnector.org/data/>. This dataset is an excerpt from the [Gapminder](http://www.gapminder.org/) data, that's [already been cleaned up to a degree](https://github.com/jennybc/gapminder). This particular dataset has 1704 observations on six variables:

* `country` a categorical variable (aka "factor") 142 levels
* `continent`, a categorical variable with 5 levels
* `year`: going from 1952 to 2007 in increments of 5 years
* `pop`: population
* `gdpPercap`: GDP per capita
* `lifeExp`: life expectancy

Let's load the data first. There are two ways to do this. You can use RStudio's menus to select a file from your computer (tools, import dataset, from text file). But that's not reproducible. The best way to do this is to save the data and the script you're using to the same place, and read in the data in using `read.csv`. It's important to tell R that the file has a header, which tells R the names of the columns. We tell this to R with the `header=TRUE` argument. 

Once we've loaded it we can type the name of the object itself (`gm`) to view the entire data frame. *Note: doing this with large data frames can cause you trouble.*


```r
gm <- read.csv("data/gapminder.csv", header=TRUE)

# Alternatively, read directly from the web:
# gm <- read.csv(url("http://bioconnector.org/data/gapminder.csv"), header=TRUE)

class(gm)
gm
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


```r
View(gm)
```

## The dplyr package

The [dplyr package](https://github.com/hadley/dplyr) is a relatively new R package that makes data manipulation fast and easy. It imports functionality from another package called magrittr that allows you to chain commands together into a pipeline that will completely change the way you write R code such that you're writing code the way you're thinking about the problem.

### tbl_df

Let's go back to our gapminder data. It's currently stored in an object called `gm`. Remember that `gm` is a `data.frame` object if we look at its `class()`. Also, remember what happens if we try to display `gm` on the screen? It's too big to show us and it fills up our console.


```r
class(gm)
gm
```

The first really nice thing about the dplyr package is the `tbl_df` class. A `tbl_df` is basically an improved version of the `data.frame` object. The main advantage to using a `tbl_df` over a regular data frame is the printing: `tbl` objects only print a few rows and all the columns that fit on one screen, describing the rest of it as text. We can turn our gm data frame into a `tbl_df` using the `tbl_df()` function. Let's do that, and reassign the result back to gm. Now, if we take a look at gm's class, we'll see that it's still a data frame, but it's also now a tbl_df. If we now type the name of the object, it will by default only print out a few lines. If this was a "wide" dataset with many columns, it would also not try to show us everything.


```r
library(dplyr)
gm <- tbl_df(gm)
class(gm)
gm
```

You don't have to turn all your data.frame objects into tbl_df objects, but it does make working with large datasets a bit easier.

---

**EXERCISE**

xx download new gapminder data, BMI? turn it into a data frame, and store it as a tbl_df. How many rows? Columns? Save it because we're going to use it later for joins.

---

### dplyr verbs

The dplyr package gives you a handful of useful **verbs** for managing data. On their own they don't do anything that base R can't do. 

1. filter
2. select
3. mutate
4. summarize
5. group_by
6. _more later..._

They all take a data.frame or tbl_df as their input for the first argument.

#### filter()

First, this is from the introduction lecture, if you want to filter **rows** of the data where some condition is true, use the `filter()` function. The first argument is the data frame or tbl, and the second argument is the condition you want to satisfy.


```r
# Show only stats for the year 1982
filter(gm, year==1982)

# Show only stats for the US
filter(gm, country=="United States")

# Show only stats for American contries in 1997
filter(gm, continent=="Americas" & year==1997)

# Show only stats for countries with per-capita GDP of less than 300 OR a 
# life expectancy of less than 30. What happened those years? 
filter(gm, gdpPercap<300 | lifeExp<30)
```

#### select()

The `filter()` function allows you to return only certain rows matching a condition. The `select()` function lets you subset the data and restrict to a number of columns. The first argument is the data, and subsequent arguments are the columns you want. Let's just get the year and the population variables.


```r
select(gm, year, pop)
```

#### mutate()

The `mutate()` function adds new columns to the data. Remember, the variable in our dataset is GDP per capita, which is the total GDP divided by the population size for that country, for that year. Let's mutate this dataset and add a column called gdp:


```r
mutate(gm, gdp=pop*gdpPercap)
```

Mutate has a nice little feature too in that it's "lazy." You can mutate and add one variable, then continue mutating to add more variables based on that variable. Let's make another column that's GDP in billions.


```r
mutate(gm, gdp=pop*gdpPercap, gdpBil=gdp/1e9)
```

#### summarize()

The `summarize()` function summarizes multiple values to a single value. On its own the `summarize()` function doesn't seem to be all that useful.


```r
summarize(gm, mean(pop))
summarize(gm, meanpop=mean(pop))
summarize(gm, n())
summarize(gm, n_distinct(country))
```

#### group_by()

We saw that `summarize()` isn't that useful on its own. Neither is `group_by()` All this does is takes an existing tbl and coverts it into a grouped tbl where operations are performed by group.


```r
gm
group_by(gm, continent)
class(group_by(gm, continent))
```

The real power comes in where `group_by()` and `summarize()` are used together. Let's take the same grouped tbl from last time, and pass all that as an input to summarize, where we get the mean population size. We can also group by more than one variable.


```r
summarize(group_by(gm, continent), mean(pop))

group_by(gm, continent, year)
summarize(group_by(gm, continent, year), mean(pop))
```

## The almighty pipe: **%>%**

This is where things get awesome. The dplyr package imports functionality from the [magrittr](https://github.com/smbache/magrittr) package that lets you _pipe_ the output of one function to the input of another, so you can avoid nesting functions. It looks like this: `%>%`. You don't have to load the magrittr package to use it since dplyr imports its functionality when you load the dplyr package.

Here's the simplest way to use it. Think of the `head()` function. It expects a data frame as input, and the next argument is the number of lines to print. These two commands are identical:


```r
head(gm, 5)
```

```
## Source: local data frame [5 x 6]
## 
##       country continent year lifeExp      pop gdpPercap
## 1 Afghanistan      Asia 1952  28.801  8425333  779.4453
## 2 Afghanistan      Asia 1957  30.332  9240934  820.8530
## 3 Afghanistan      Asia 1962  31.997 10267083  853.1007
## 4 Afghanistan      Asia 1967  34.020 11537966  836.1971
## 5 Afghanistan      Asia 1972  36.088 13079460  739.9811
```

```r
gm %>% head(5)
```

```
## Source: local data frame [5 x 6]
## 
##       country continent year lifeExp      pop gdpPercap
## 1 Afghanistan      Asia 1952  28.801  8425333  779.4453
## 2 Afghanistan      Asia 1957  30.332  9240934  820.8530
## 3 Afghanistan      Asia 1962  31.997 10267083  853.1007
## 4 Afghanistan      Asia 1967  34.020 11537966  836.1971
## 5 Afghanistan      Asia 1972  36.088 13079460  739.9811
```

So what? 

Now, think abou this about this for a minute. What if we wanted to get the life expectancy and GDP averaged across all Asian countries for each year? (See top of image). Mentally we would do something like this:

0. Take the `gm` dataset
0. `mutate()` it to add raw GDP
0. `filter()` it to restrict to Asian countries only
0. `group_by()` year
0. and `summarize()` it to get the mean life expectancy and mean GDP.

But in code, it gets ugly. First, `mutate` the data to add GDP.


```r
mutate(gm, gdp=gdpPercap*pop)
```

Wrap that whole command with `filter()`.


```r
filter(mutate(gm, gdp=gdpPercap*pop), continent=="Asia")
```

Wrap that again with `group_by()`:


```r
group_by(filter(mutate(gm, gdp=gdpPercap*pop), continent=="Asia"), year)
```

Finally, wrap everything with `summarize()`:


```r
summarize(
  group_by(
    filter(
      mutate(gm, gdp=gdpPercap*pop), 
    continent=="Asia"), 
  year), 
mean(lifeExp), mean(gdp))
```

Now compare that with the mental process of what you're actually trying to accomplish. The way you would do this without pipes is completely inside-out and backwards from the way you express in words and in thought what you want to do. The pipe operator `%>%` allows you to pass data frame or tbl objects from one function to another, so long as those functions expect data frames or tables as input.

<img src="{{site.baseurl}}/assets/gmdplyr.png", width=750>

This is how we would do that in code. It's as simple as replacing the word "then" in words to the symbol `%>%` in code.


```r
gm %>%
  mutate(gdp=gdpPercap*pop) %>%
  filter(continent=="Asia") %>%
  group_by(year) %>%
  summarize(mean(lifeExp), mean(gdp))
```

