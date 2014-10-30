---
layout: page
---



# Introduction to R

*[back to course contents](..)*

The first part of this workshop will demonstrate very basic functionality in R, including functions, functions, vectors, creating variables, getting help, subsetting, data frames, plotting, and reading/writing files.

[Link to slides](https://speakerdeck.com/stephenturner/introduction-to-r-for-life-scientists).

[Link to R Cheat Sheet](http://dx.doi.org/10.6084/m9.figshare.1080756).

## Before coming

{% include setup-r.md %}

## RStudio

Let's start by learning about RStudio. **R** is the underlying statistical computing environment, but using R alone is no fun. **RStudio** is a graphical integrated development environment that makes using R much easier.

- Panes in RStudio: I set up my window to have the editor in the top left, console top right, environment/history on the bottom left, and plots/help on the bottom right. 
- Code that you type into the console is code that R executes. From here forward we will use the editor window to write a script that we can save to a file and run it again whenever we want to. We usually give it a `.R` extension, but it's just a plain text file. If you want to send commands from your editor to the console, use `CMD`+`Enter` (`Ctrl`+`Enter` on Windows).
- Anything after a `#` sign is a comment. Use them liberally to *comment your code*.

## Basic operations

R can be used as a glorified calculator. Try typing this in directly into the console. Make sure you're typing into into the editor, not the console, and save your script. Use the run button, or press `CMD`+`Enter` (`Ctrl`+`Enter` on Windows).


```r
2 + 2
5 * 4
2^3
```

R Knows order of operations and scientific notation.


```r
2 + 3 * 4/(5 + 3) * 15/2^2 + 3 * 4^2
50000
```

However, to do useful and interesting things, we need to assign *values* to *objects*. To create objects, we need to give it a name followed by the assignment operator `<-` and the value we want to give it:


```r
weight_kg <- 55
```

`<-` is the assignment operator. Assigns values on the right to objects on the left, it is like an arrow that points from the value to the object. Mostly similar to `=` but not always. Learn to use `<-` as it is good programming practice. Using `=` in place of `<-` can lead to issues down the line.

Objects can be given any name such as `x`, `current_temperature`, or `subject_id`. You want your object names to be explicit and not too long. They cannot start with a number (`2x` is not valid but `x2` is). R is case sensitive (e.g., `weight_kg` is different from `Weight_kg`). There are some names that cannot be used because they represent the names of fundamental functions in R (e.g., `if`, `else`, `for`, see [here](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Reserved.html) for a complete list). In general, even if it's allowed, it's best to not use other function names, which we'll get into shortly (e.g., `c`, `T`, `mean`, `data`, `df`, `weights`). In doubt check the help to see if the name is already in use. It's also best to avoid dots (`.`) within a variable name as in `my.dataset`. It is also recommended to use nouns for variable names, and verbs for function names.

When assigning a value to an object, R does not print anything. You can force to print the value by typing the name:


```r
weight_kg
```

Now that R has `weight_kg` in memory, we can do arithmetic with it. For instance, we may want to convert this weight in pounds (weight in pounds is 2.2 times the weight in kg):


```r
2.2 * weight_kg
```

We can also change a variable's value by assigning it a new one:


```r
weight_kg <- 57.5
2.2 * weight_kg
```

This means that assigning a value to one variable does not change the values of other variables. For example, let's store the animal's weight in pounds in a variable.


```r
weight_lb <- 2.2 * weight_kg
```

and then change `weight_kg` to 100.


```r
weight_kg <- 100
```

What do you think is the current content of the object `weight_lb`? 126.5 or 220?

You can see what objects (variables) are stored by viewing the Environment tab in Rstudio. You can also use the `ls()` function. You can remove objects (variables) with the `rm()` function. You can do this one at a time or remove several objects at once.


```r
ls()
rm(weight_lb, weight_kg)
ls()
weight_lb  # oops! you should get an error because weight_lb no longer exists!
```

---

**EXERCISE**

What are the values after each statement in the following?


```r
mass <- 50  # mass?
age <- 30  # age?
mass <- mass * 2  # mass?
age <- age - 10  # age?
mass_index <- mass/age  # massIndex?
```

---

## Functions

R has built-in functions.


```r
# Notice that this is a comment.  Anything behind a # is 'commented out' and
# is not run.
sqrt(144)
log(1000)
```

Get help by typing a question mark in front of the function's name, or `help(functionname)`:

```
help(log)
?log
```

Note syntax highlighting when typing this into the editor. Also note how we pass *arguments* to functions. The `base=` part inside the parentheses is called an argument, and most functions use arguments. Arguments modify the behavior of the function. Functions some input (e.g., some data, an object) and other options to change what the function will return, or how to treat the data provided. Finally, see how you can *next* one function inside of another (here taking the square root of the log-base-10 of 1000).


```r
log(1000)
log(1000, base = 10)
log(1000, 10)
sqrt(log(1000, base = 10))
```

---

**EXERCISE**

See `?abs` and calculate the square root of the log-base-10 of the absolute value of `-4*(2550-50)`. Answer should be `2`.

---

## Vectors and classes

Let's create some numeric vectors. Vectors (aka "arrays" in Perl, "lists" in Python) are single *objects* containing an ordered collection of *elements*. A simple vector is a numeric vector, a single *object* containing several numbers. Here let's display a few vectors. We can also do vector arithmetic. When printing vectors to the screen that have lots of elements, notice that the bracketed number in the gutter of the output is just a counter indexing the number of elements in the vector.


```r
# Some simple numeric vectors:
1:5
6:10
1:5 + 6:10
1:100
```

We can also create arbitrary vectors with the `c()` function (short for "combine").


```r
c(1, 2, 5)
c(1:5, 11:15)
```

What if we wanted to create a vector from 2 to 10 by 2's? What about 2 to 200 by 4's? This might be useful for setting up an experiment where every other sample is an experimental group and every other is a control.


```r
c(2, 4, 6, 8, 10)

# Get some help with the seq() function, then create a vector from 2 to 200
# by 2s.  Notice how the seq() function works -- the `to` argument will
# never be exceeded.
help(seq)
seq(from = 2, to = 200, by = 4)
```

You can assign this vector of values to an object, just like you would for one item. For example we can create a vector of animal weights:


```r
weights <- c(50, 60, 65)
weights
```

A vector can also contain characters:


```r
animals <- c("mouse", "rat", "dog")
animals
```

There are many functions that allow you to inspect the content of a vector. `length()` tells you how many elements are in a particular vector:


```r
length(weights)
length(animals)
```

`class()` indicates the class (the type of element) of an object:


```r
class(weights)
class(animals)
```

---

**EXERCISE**

- Use the `c()` function to create/assign a new object that combines the `weights` and `animals` vectors into a single vector called `combined`.
- What happened to the numeric values? *Hint*: What's the `class()` of `combined`?
- Why do you think this happens?

---

The function `str()` provides an overview of the structure of an object and the elements it contains. It is a really useful function when working with large and complex objects:


```r
str(weights)
str(animals)
```

You can add elements to your vector simply by using the `c()` function:


```r
weights
weights <- c(weights, 90)  # adding at the end
weights <- c(30, weights)  # adding at the beginning
weights
```

What happens here is that we take the original vector `weights`, and we are adding another item first to the end of the other ones, and then another item at the beginning. We can do this over and over again to build a vector or a dataset. When you're programming this may be useful to autoupdate results that we are collecting or calculating.

Certain *functions* operate only on certain *classes* of object. Here, `weights` is a `numeric` vector. The built-in `sum()` function will operate on numeric objects, but not characters.


```r
sum(weights)
sum(animals)
```

---

**EXERCISE**

Sum the integers 1 through 100 and 501 through 600 (e.g. 1+2+...+99+100+501+502+...+599+600)

---

## Slicing/indexing vectors

Let's create a vector of 50 integers going from 101 to 150. We can access certain elements of that vector by putting the element's *index(es)* in square brackets. E.g., `x[1]` will return the first element in vector `x`. Calling `x[c(3,5)]` will access the third and fifth elements. Calling `x[1:10]` will return the first ten elements of `x`.

*Special note: R indexes vectors starting at 1. This is different from many other languages, including Perl and Python, which index starting from zero.*


```r
# Create the vector.
x <- 101:150

# Get the first element.
x[1]

# Get the 42nd element.
x[42]

# Get the 20th through the 25th elements.
x[20:25]

# If you try to access elements that don't exist, you'll return missing
# values.  Missing values are represented as NA
x[45:55]  #NA is missing value!
```

## Data Frames

We use **data frames** to store heterogeneous tabular data in R: tabular, meaning that individuals or observations are typically represented in rows, while variables or features are represented as columns; heterogeneous, meaning that columns/features/variables can be different classes (on variable, e.g. age, can be numeric, while another, e.g., cause of death, can be text).

Later on we'll read in data from a text file into a data frame object using one of the functions `read.table()` for text files or `read.csv()` for comma-separated tables. But for now, let's use a built-in data frame called `mtcars`. This data was extracted from the 1974 Motor Trend US magazine, and comprises fuel consumption and 10 aspects of automobile design and performance for 32 automobiles (1973–74 models). We can load this built-in data with `data(mtcars)`. By the way, running `data()` without any arguments will list all the available built-in datasets included with R.

Let's load the data first. Type the name of the object itself (`mtcars`) to view the entire data frame. *Note: doing this with large data frames can cause you trouble.*


```r
data(mtcars)
class(mtcars)
mtcars
```

### Inspecting data.frame objects

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
    * `rownames()`: returns the row names
    * `str()`: structure of the object and information about the class, length and content of each column
    * `summary()`: works differently depending on what kind of object you pass to it. Passing a data frame to the `summary()` function prints out some summary statistics about each column (min, max, median, mean, etc.)


```r
head(mtcars)
tail(mtcars)
dim(mtcars)
nrow(mtcars)
ncol(mtcars)
colnames(mtcars)
rownames(mtcars)
str(mtcars)
summary(mtcars)
```

### Accessing variables & subsetting data frames

We can access individual variables within a data frame using the `$` operator, e.g., `mydataframe$specificVariable`. Let's print out the number of cylinders for every car, and calculate the average miles per gallon for ever car in the dataset (using the built-in `mean()` function).


```r
# display the number of cylinders for each car.
mtcars$cyl
# first display MPG for all vehicles, then calculate the average.
mtcars$mpg
mean(mtcars$mpg)
```

We can also use the `subset()` function to return a subset of the data frame that meets a specific condition. 

1. The first argument is the data frame you want to subset, e.g. `subset(mtcars...`.
2. The second argument is a condition you must satisfy, e.g. `subset(mtcars, cyl==6)`. If you want to satisfy *all* of multiple conditions, you can use the "and" operator, `&`. The "or" operator `|` (the pipe character, usually shift-backslash) will return a subset that meet *any* of the conditions.
    * `==`: Equal to
    * `!=`: Not equal to
    * `>`, `>=`: Greater than, greater than or equal to
    * `<`, `<=`: Less than, less than or equal to

Try it out:


```r
# Return only cars with 6 cylinder engines.
subset(mtcars, cyl == 6)
# Return only cars having more than 6 cylinders **and** an engine
# displacement volume less than 300cc.
subset(mtcars, cyl > 6 & disp < 300)
# Return only the cars that get at least 20 miles per gallon or have a
# displacement volume of less than 100cc.
subset(mtcars, mpg >= 20 | disp < 100)
```

Finally, take a look at the class of what's returned by a `subset()` function. The `subset()` function takes a data.frame and returns a data.frame. You can operate on this new data.frame just as you would any other data.frame using the `$` operator. E.g., print the MPG for all the 6 cylinder vehicles:


```r
subset(mtcars, cyl == 6)$mpg
```

---

**EXERCISE**

1. Print out the dataset cars that have greater than or equal to 6 cylinders *and* get at least 15 miles per gallon.
2. What's the mean engine displacement of these vehicles?



---

### with()

The `with()` function is particularly helpful. Let's say you wanted to compute some (senseless) value by computing the MPG times the number of cylinders divided by the car's displacement. You could access the dataset's variables using the `$` notation, or you could use `with()` to temporarily *attach* the data frame, and call the variables directly, as if they were just vectors hanging out in your workspace. The first argument to `with()` is the name of the data frame, and the second argument is all the stuff you'd like to do with the particular features in that data frame.

Try typing the following commands:


```r
# Display the number of cylinders.
mtcars$cyl
with(mtcars, cyl)

# Compute the senseless value described above. Both return the same results.
mtcars$mpg * mtcars$cyl/mtcars$disp
with(mtcars, mpg * cyl/disp)
```

---

**EXERCISE**

Using the `with()` and `subset()` functions, compute the average value of the ratio of the engine displacement (`disp`) divided by the engine's horsepower (`hp`) for all the 8-cylinder vehicles.


```
## [1] "Answer: 1.754801425962"
```

---

## Plotting basics

Plotting a single numeric variable goes down the rows and plots a value on the y-axis for each observation (index) in the data frame.


```r
plot(mtcars$mpg)
```

This isn't a very useful figure. More appropriate might be a histogram. We can try to let R decide how many breaks to insert in the histogram, or we can set that manually. We can also set the color of the bars.


```r
hist(mtcars$mpg)
hist(mtcars$mpg, breaks = 10)
hist(mtcars$mpg, breaks = 10, col = "black")
```

We can create a scatterplot between two variables with `plot(varX, varY)`.


```r
# This would also work, but let's use with().  plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))
```

There are hundreds of plotting parameters you can use to make your plot look exactly like you want. Let's use a solid-filled point instead of an open circle with the `pch=` argument (point character), color the points red with the `col=` argument, give it a title by passing a character object to the `main=` argument, and change the x and y axis titles with the `xlab=` and `ylab=` arguments, respectively. Let's go through this one step at a time.


```r
with(mtcars, plot(disp, mpg, pch = 16))
with(mtcars, plot(disp, mpg, pch = 16, col = "red"))
with(mtcars, plot(disp, mpg, pch = 16, col = "red", main = "MPG vs Displacement"))
with(mtcars, plot(disp, mpg, pch = 16, col = "red", main = "MPG vs Displacement", 
    ylab = "Fuel Economy (MPG)", xlab = "Displacement (cu. in.)"))
```

Notice how on that last line I broke the command up into two lines for better readability. I broke the command at the comma separating arguments, and indented the following line for readability.

With plotting parameters, **Google is your friend.** 

* Forget what each point character represents? Google _R pch_.
* Forget the names of R's colors? Google _R colors_. Want to learn more about color schemes in R? Google _RColorBrewer_.
* Try googling _R graphical parameters_.

---

**EXERCISE**

Plot horsepower (y-axis) vs displacement (x-axis) for vehicles with more than 4 cylinders. Give the graph a title and label the axes. Make the points solid (hint, `pch=16`) blue (hint, `col="blue"`) circles.



---

## Reading in / writing out data

First, lets create a small dataset consisting of only 8 cylinder cars.


```r
mtcars_8cyl <- subset(mtcars, cyl == 8)
mtcars_8cyl
```

Next, check what your working directory is with `getwd()` with no arguments, and look up some help for `write.table()` and `write.csv()`.


```r
getwd()
help(write.table)
help(write.csv)
```

Using RStudio, go to the Session menu, and select the directory (folder) you want to work from under the "Set Working Directory" menu. You can also do this manually with the `setwd()` command.


```r
getwd()
setwd("~/Desktop/R")
```

Once you've set your working directory either using RStudio or on the command line, save the new reduced data frame to a comma-separated file called `cars8.csv` using the `write.csv()` function.


```r
write.csv(mtcars_8cyl, file = "cars8.csv")
```

Data can be loaded using the Tools -- Import Dataset -- From text file menu in R studio. Or you can also load a dataset manually using `read.table()` or `read.csv()`. First, read the help on these functions:


```r
help(read.table)
help(read.csv)
```

Here let's remove the dataset, and re-import it into an object called cars8 from the file we just saved.


```r
rm(mtcars_8cyl)
mtcars_8cyl
cars8 <- read.table(file = "cars8.csv", header = TRUE, sep = ",", row.names = 1)
cars8
rm(cars8)
cars8 <- read.csv(file = "cars8.csv", header = TRUE, row.names = 1)
cars8
```
