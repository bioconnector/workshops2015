# Notes for teaching intro R

## Links

- Lesson source: [Rmd](r-intro.Rmd), [md](r-intro.md)
- Rendered lesson on [bioconnector.org](http://bioconnector.org/workshops/lessons/r/r-intro/)
- Exercises: [source md](r-intro-ex.md), [PDF handout](r-intro-ex.md.pdf)
- Data: [bioconnector.org/data](http://bioconnector.org/data/) (gapminder.csv)

## Teaching notes

### RStudio

- I usually start the day off with a short slideshow to motivate participants to learn R ([link to slides](https://speakerdeck.com/stephenturner/introduction-to-r-for-life-scientists)). I'm happy to do this section.
- When I first dive into hands on, I tell participants to launch RStudio. I make clear the difference between R (the underlying statistical computing engine/environment) versus RStudio (the integrated development environment). Emphasize that it's important to launch RStudio, not R.
- The first thing that I do is change my pane layout. I set up my window to have the editor in the top left, console top right, environment/history on the bottom left, and plots/help on the bottom right. Emphasize that participants _don't have to_ change their pane layout to match mine. This layout is my preference, but if they want their panes to look like mine, I highly recommend it. Do this through Tools -- Global Options -- Pane Layout.
- Project management
    - The next thing I do is tell participants to start a new project. File -- New Project -- New Directory -- Empty Project. Create the project somewhere easy to find (e.g., Desktop), and call it something easy to remember (workshop).
    - Completely quit/close RStudio at this point. 
    - Navigate to the project folder you created and show participants the .Rproj file you created. Emphasize that creating an R project allows you to open R and RStudio _running in the project folder_. This is important for reproducibility and making loading data easier - you'll always want to save your R code and your data in the same place. 
    - Emphasize again that the .Rproj (Rproject file) isn't code or data - it's only a file that launches R _working in the project directory_. 
    - Open your .Rproject file and continue.
- At this point I'll type a command or two into the console but then I immediately tell participants not to use the console. I open up a script (File -- New -- Script), type those commands into a script, and save the script. Emphasize how this script is just a text file (optionally with a .R extention, nothing special, just text). 
- Here is where I demonstrate the keyboard shortcut to execute a command (CMD-Enter on Mac, Ctrl-Enter on Windows). I never use the "Run line" button - I always use the shortcut.
- If you use the keyboard shortcut for `<-`, you'll want to explain that.
- If you talk at length about keyboard shortcuts, you might mention in the help menu there's an item that displays all the shortcuts.

### Things people get hung up on

- `<-` and `=`: both can be used as assignment operators. I almost always use `<-`, but `=` works _most_ of the time. Pick a style and stick with it.
- In function calls, why do you use `=` instead of `<-` to pass parameters to arguments?
- `==` instead of `=` for comparisons.
- Operating on a variable versus reassignment. E.g., 
  ```r
  x <- 3
  x + 5 # doesn't change x
  x <- x+5 # does change x
  ```
  This is especially confusing when we get to dplyr / chaining in later lessons.
