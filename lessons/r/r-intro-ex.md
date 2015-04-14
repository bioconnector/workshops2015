# Exercises: Introduction to R

**Exercise 1**

What are the values after each statement in the following?


```r
mass <- 50  # mass?
age <- 30  # age?
mass <- mass * 2  # mass?
age <- age - 10  # age?
mass_index <- mass/age  # massIndex?
```

----

**Exercise 2**

See `?abs` and calculate the square root of the log-base-10 of the absolute value of `-4*(2550-50)`. Answer should be `2`.

----

**Exercise 3**

- Use the `c()` function to create/assign a new object that combines the `weights` and `animals` vectors into a single vector called `combined`.
- What happened to the numeric values? *Hint*: What's the `class()` of `combined`?
- Why do you think this happens?

----

**Exercise 4**

Sum the integers 1 through 100 and 501 through 600 (e.g. 1+2+...+99+100+501+502+...+599+600)

----

**Exercise 5**

1. What country and what years had a low GDP (<500) but high life expectancy (>50)?
2. What's the average GDP for Asian countries in 2002? How does that compare to European countries in the same year? To the Americas?


----

**Exercise 6**

Using the `with()`, do the following: 

1. Compute the average GDP in billions for all Asian countries in 2007. 
2. Do the same for Europe in 2007.

_Hint:_ GDP per capita is the GDP divided by the population size. So to get GDP, you'd multiple `gdpPercap*pop`. To get that in billions, divide by 1,000,000, or more easily expressed in R using scientific notation: `1e9`.

----

**Exercise 7**

Plot GDP in trillions (`gdpPercap*pop/1e9`) on the y-axis versus population size in millions on the x-axis for all countries in the Americas. Use solid (`pch=16`) "blue" points, and give the plot a title and legends.
