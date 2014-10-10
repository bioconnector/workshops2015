---
layout: page
---

# Introduction to the UNIX shell

Modified from Tracy Teal's lessons for Data Carpentry with original contributions from Paul Wilson, Milad Fatenejad, Sasha Wood and Radhika Khetani for Software Carpentry.

Learn more:

- <http://datacarpentry.org/>
- <http://software-carpentry.org/>

## Objectives

-	What is the shell?
-	How do you access it?
- Connecting to cloud computing resources
-	How do you use it?
  -	Getting around the Unix file system
  -	looking at files
  -	manipulating files
  -	automating tasks
-	What is it good for?
-	Where are resources where I can learn more? (because the shell is awesome)

## What is the shell?

The *shell* is a program that presents a command line interface which allows you to control your computer using commands entered with a keyboard instead of controlling graphical user interfaces (GUIs) with a mouse/keyboard combination.

*Why should you care?*

- For 99% of bioinformatics tools, you have to use the shell (the command line). There is no graphical interface.
- The shell gives you *power*. The command line gives you the power to do your work more efficiently and more quickly. When you need to do things tens to hundreds of times, knowing how to use the shell is transformative.
- You have to use the shell to connect to remote computers.

### Automation

Have 10,000,000 files to rename, read in, analyze, and visualize? It's easy to automate things with the shell.

![Automation](../img/geek_vs_nongeek.png)


### How to access the shell

You access the shell through a program called a Terminal. We're going to use a terminal to connect to the shell of a remote computer. The terminal program is built-in on Mac and Linux. For Windows, you'll have to download a separate program called a terminal emulator that will allow you to connect to remote computers.

#### Mac

On Mac the shell is available through Terminal:

**Applications -> Utilities -> Terminal**

Go ahead and drag the Terminal application to your Dock for easy access.

#### Windows

On Windows machines we'll be using a terminal emulator called [PuTTY](http://the.earth.li/~sgtatham/putty/latest/x86/putty.exe).


### More resources on the shell

Cheat sheets:

- <http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/>

Web sites where you can see what the different components of a shell command are doing:

- [explainshell.com](http://explainshell.com)
- [commandlinefu.com](http://www.commandlinefu.com)


## Connecting to cloud computing resources

Cloud computing means different things to different people. There's Software-as-a-Service (Saas) like Gmail, Facebook, etc., that most people think of when they think of cloud computing. Here, we're strictly talking about Infrastructure-as-a-Service (IaaS). That is, we're using (or *renting*) computing infrastructure that we don't *own*.

fixme:
- fixme: add content about amazon on separate 00-lesson?
- fixme: launching instance w/o keys
- fixme: connecting to instance w/ username/password (terminal on mac, putty on windows)
- fixme: file transfer on mac/windows with cyberduck

## Shell basics

### Get some data files to work with

We will spend most of our time learning about the basics of the shell by manipulating some experimental data. Now we're going to download the data for the tutorial. For this you'll need internet access, because you're going to get it off the web.

Open the shell and type the command:

```
git clone https://github.com/bioconnector/workshops.git
```

This command will grab all of the data needed for this workshop. It's using something called git that's used for version control, but we won't talk about that here.

### Moving around and listing files

We'll start out by moving around the file system and listing files. Today we're going to go through using the command line. These commands are in the README.md file and on your handout (fixme).

Let's go in to that directory we just downloaded:

```
cd workshops
```

`cd` stands for 'change directory'. A *directory* is the same thing as a *folder*. We just call folders *directories* in the Linux/UNIX world. Here we just entered the workshops directory. If we were doing this on a graphical shell we would have double-clicked on a little folder icon. It's the same idea.

In this directory, there should be some things we just downloaded. Let's check. Type:

```
ls
```

`ls` stands for 'list' and it lists the contents of a directory. This woudld be what you would see in the folder if you were doing this graphically on a desktop.

Blue things are directories, white things are files.

Now, let's go look in the 'data' directory in the 'shell' lesson. It's nested a few directories deep. To get to it, let's enter and list each directory like so:

```
cd lessons
ls
cd shell
ls
cd data
ls
```

In there, all mixed up together are regular files, directories, and an executable program. If we want to know which is which, we can type:

```
ls -F
```

Things with a `/` after it is a directory.  
Things with a `*` after them are programs.  
It there's nothing there it's a regular file.

You can also use the command:

```
ls -l
```

to see whether items in a directory are files or directories. It gives a lot more information too, such as the size of the file, who owns the file, etc.

So, we can see that we have several files, directories and a program. Great!

### Options

Most programs take additional options that control their exact behavior. For example, `-F` and `-l` are options for `ls`. The `ls` program, like many programs, take a lot of options. But how do we know what the options are to particular commands?

Most commonly used shell programs have a manual. Let's open the manual page for `ls`.

You can access the manual using the `man` program.

```
man ls
```

Space key goes forward  
Or use the arrow keys to scroll up and down.  
When you are done reading, just hit `q` to quit.

Programs that are run from the shell can get extremely complicated. To see an example, open up the manual page for the `find` program. No one can possibly learn all of these arguments, of course. So you will probably find yourself referring back to the manual page frequently.

### The Unix directory file structure (a.k.a. where am I?)

As you've already just seen, you can move around in different directories or folders at the command line.

When you're working with bioinformatics programs, you're working with your data and it's key to be able to have that data in the right place and make sure the program has access to the data. Many of the problems people run in to with command line bioinformatics programs is not having the data in the place the program expects it to be.

Let's draw out what we just did, and some of the other files and folders we could have clicked on.

This is called a hierarchical file system structure, like an upside down tree with root (/) at the base that looks like this.

![Unix](../img/unix_filesystem.png)

That (/) at the base is often also called the 'top' level.

When you are working at your computer or log in to a remote computer, you are on one of the branches of that tree, your home directory: **/home/username**.

If we type `cd` by itself:

```
cd
```

This puts you in your home directory. That's **/home/username**

---

**EXERCISE

-	Using `cd` and `ls`, go in to the 'workshops/shell/data' directory and list its contents.
-	How many files, how many directories and how many programs are there?

---

### Where am I?

Let's also check to see where we are. Sometimes when we're wandering around in the file system, it's easy to lose track of where we are and get lost.

If you want to know what directory you're currently in, type:

```
pwd
```

This stands for 'print working directory'. That's the directory you're currently working in, and it's "printed" to the screen.

What if we want to move back up and out of the `data` directory? To go 'back up a level' we need to use `..`

```
cd ..
```

Now do `ls` and `pwd`. See now that we went back up in to the 'shell' directory. `..` means "the directory above," or "the parent directory."

---

**EXERCISE**

Let's go on a file hunt. Move around in the "shell/data/hidden" directory and try to find the file "youfoundit.txt."

---

### Examining the contents of other directories

By default, the `ls` commands lists the contents of the working directory (i.e. the directory you are in). You can always find the directory you are in using the `pwd` command. However, you can also give `ls` the names of other directories to view. Navigate to the home directory if you are not already there using `cd` by itself, then `ls` the contents of the "workshops/lessons/shell" directory:

```
cd
ls workshops/lessons/shell
```

This listed the contents of workshops/lessons/shell without navigating there.

The `cd` command works the same way. Try entering:


```
cd
ls workshops/lessons/shell/data/hidden
```

and you will jump directly to `hidden` without having to go through the intermediate directories.

---

**EXERCISE**

Try finding the 'anotherfile.txt' file without changing directories.

---


### HUGE Shortcut: Tab Completion

Navigate to the home directory. Typing out directory names can waste a lot of time. When you start typing out the name of a file or directory, then hit the tab key, the UNIX shell will try to fill in the rest of the directory name. For example, enter:

```
cd
cd w<tab>
```

The shell will fill in the rest of the directory name for "workshops". Now go to `workshops/lessons/rnaseq-1day/data`

```
cd
cd wo<tab>le<tab>rna<tab>da<tab>
```

Now type `ls F3D` and hit tab twice.

```
ls uvb<tab><tab>
```

When you hit the first tab, nothing happens. The reason is that there are multiple directories in the home directory which start with `uvb`. Thus, the shell does not know which one to fill in. When you hit tab again, the shell will list the possible choices.

Tab completion can also fill in the names of programs. For example, enter `e<tab><tab>`. You will see the name of every program that starts with an `e`. One of those is `echo`. If you enter `ec<tab>` you will see that tab completion works.

### Full vs. Relative Paths

The `cd` command takes an argument which is the directory name. Directories can be specified using either a *relative path* or a *full path*. The directories on the computer are arranged into a hierarchy. The full path tells you where a directory is in that hierarchy. Navigate to the home directory. Now, enter the `pwd` command and you should see:

```
/home/username
```

which is the full name of your home directory. This tells you that you are in a directory called `username`, which sits inside a directory called `home` which sits inside the very top directory in the hierarchy. The very top of the hierarchy is a directory called `/` which is usually referred to as the *root directory*. So, to summarize: `username` is a directory in `home` which is a directory in `/`.

Now enter the following command:

```
cd /home/ubuntu/workshops/lessons/shell/data/hidden
```

This jumps to `hidden`. Now go back to the home directory (`cd`). We saw earlier that the command:

```
cd workshops/lessons/shell/data/hidden
```

had the same effect - it took us to the `hidden` directory. But, instead of specifying the full path, we specified a *relative path*. In other words, we specified the path relative to our current directory. A full path always starts with a `/`. A relative path does not.

A relative path is like getting directions from someone on the street. They tell you to "go right at the Stop sign, and then turn left on Main Street". That works great if you're standing there together, but not so well if you're trying to tell someone how to get there from another country. A full path is like GPS coordinates. It tells you exactly where something is no matter where you are right now.

You can usually use either a full path or a relative path depending on what is most convenient. If we are in the home directory, it is more convenient to just enter the relative path since it involves less typing.

Over time, it will become easier for you to keep a mental note of the structure of the directories that you are using and how to quickly navigate them.

---

**EXERCISE**

-	List the contents of the /bin directory. Do you see anything familiar in there?

---

### Saving time with shortcuts and wild cards

#### Shortcuts

There are some shortcuts which you should know about. Dealing with the home directory is very common. In the shell the tilde character, `~`, is a shortcut for your home directory.

```
ls ~
```

Try it even when you're not in your home directory:

```
cd
cd workshops
```

Then enter the command:

```
ls ~
```

This prints the contents of your home directory, without you having to type the full path. The shortcut `..` always refers to the directory above your current directory.

```
ls ..
```

prints the contents of the `/home/username/`. You can chain these together, so:

```
ls ../../
```

prints the contents of `/home/`. Finally, the special directory `.` always refers to your current directory. So, `ls`, `ls .`, and `ls ././././.` all do the same thing, they print the contents of the current directory. This may seem like a useless shortcut right now, but we'll see when it is needed in a little while.

To summarize, while you are in the `workshops` directory, the commands `ls ~`, `ls ../`, and `ls /home/ubuntu` all do exactly the same thing.

#### Our data set: FASTQ files

In this example we're going to use some RNA-seq data that we'll actually analyze later on. The data come from RNA-seq done on skin cells treated with ultraviolet light versus controls ([paper](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0093338); [data](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54413)).

We did sequencing, and what we have here are just reads that came from a 100 MB region on chromosome 4. We want to be able to look at these files and do some things with them.

#### Wild cards

Navigate to the `workshops/lessons/rnaseq-1day/data` directory from your home. This directory contains our FASTQ files and some other ones we'll need for analyses. If we type `ls`, we will see that there are several files in there. Some of the end with .fastq.gz. These are compressed fastq files.

The `*` character is a shortcut for "everything". Thus, if you enter `ls *`, you will see all of the contents of a given directory.

```
ls *fastq.gz
```

This lists every file that ends with a `fastq.gz`. This command:

```
ls /usr/bin/*.sh
```

Lists every file in `/usr/bin` that ends in the characters `.sh`.

If we wanted to list just the controls or just the uvb-treated samples, we could do this:

```
ls ctl*
ls uvb*
```

So how does this actually work? Well...when the shell sees a word that contains the `*` character, it automatically looks for filenames that match the given pattern. In this case, it identified 3 such files. Then, it replaced the `ctl*` with the list of files, separated by spaces.

Because we'll use these files later, let's go ahead and extract or uncompress these files. To compres a text file we would use `gzip`. To extract, we'll use `gunzip`. Let's extract all the files that end in `.gz`. This will extract both the fastq files and the counts.txt file as well:

```
gunzip *.gz
```

Now, time for some more practice with wildcards.

---

**EXERCISE**

Do each of the following using a single `ls` command without navigating to a different directory.
  -	List all of the files in `/bin` that start with the letter 'c
  -	List all of the files in `/bin` that contain the letter 'a'
  -	List all of the files in `/bin` that end with the letter 'o'

---

### Command History

You can easily access previous commands. Hit the up arrow. Hit it again. You can step backwards through your command history. The down arrow takes your forwards in the command history.

Control-C will cancel the command you are writing, and give you a fresh prompt.

You can also review your recent commands with the `history` command to see a numbered list of commands you've run.

```
history
```

You can reuse one of these commands directly by referring to the number of that command. If your history looked like this:

```
259  ls *
260  ls /usr/bin/*.sh
261  ls ctl*
```


then you could repeat command `#260` by simply entering:

```
!260
```

---

**EXERCISE**

Find the line number in your history for the last exercise (listing files in /bin) and reissue that command.

---

### Examining Files

We now know how to switch directories, run programs, and look at the contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to print out all of the contents using the program `cat`. Enter the following command:

```
cat ctl1.fastq
```

This prints out the contents of the `ctl1.fastq` file.


---

**EXERCISE**

-	Print out the contents of the `workshops/lessons/rnaseq-1day/data/coldata.csv` file. What does this file contain?
-	Without changing directories, use one short command to print the contents of all of the files in the `~/workshops/lessons/shell/` directory.

---

Make sure we're in the right place for the next set of the lessons. We want to be in the rnaseq data directory. Check if you're there with `pwd` and if not navigate there. One way to do that would be

```
cd ~/workshops/lessons/rnaseq-1day/data
```

`cat` is a terrific program, but when the file is really big, it can be annoying to use.

The program, `less`, is useful when files are big and you want to be able to scroll through them.

```
less F3D0_S188_L001_R1_001.fastq
```

`less` opens the file, and lets you navigate through it. The commands
are identical to the `man` program.


**Some commands in `less`**

| Key   | Action              |
|-------|---------------------|
| space | go forward          |
| b     | go backwarsd        |
| g     | go to the beginning |
| G     | go to the end       |
| q     | quit                |

`less` also gives you a way of searching through files. Just hit the "/" key to begin a search. Enter the name of the word you would like to search for and hit enter. It will jump to the next location where that word is found. Try searching `ctl1.fasta` for the the sequence "GATTACA". If you hit "/" then "enter", `less` will just repeat the previous search. `less` searches from the current location and works its way forward. If you are at the end of the file and search for "GATTACA", `less` will not find it. You need to go to the beginning of the file with the `g` key in less, and start the search from there.

Remember, the `man` program actually uses `less` internally and therefore uses the same commands, so you can search documentation using "/" as well!

There's another way that we can look at files, and in this case, just look at part of them. This can be particularly useful if we just want to see the beginning or end of the file, or see how it's formatted.

The commands `head` and `tail` let you look at the beginning and end of a file respectively.

```
head ctl1.fastq
tail ctl1.fastq
```

The `-n` option to either of these commands can be used to print the first or last `n` lines of a file. If we want to see the first read in the file, try this:

```
head -n 4 ctl1.fastq
tail -n 4 ctl1.fastq
```

### Searching files

We showed a little how to search within a file using `less`.

We can search within files without even opening them, using `grep`. Grep is a command-line utility for searching plain-text data sets for lines matching a string or regular expression.

Let's find all the lines in ctl1.fastq that contain the sequence motif "GATTACA."

```
grep GATTACA ctl1.fastq
```

We get back just the sequence line, but what if we wanted all four lines, the whole part of that FASTQ sequence, back instead.

```
grep -B 1 -A 2 GATTACA ctl1.fastq
```

The `-A` flag stands for "after match" so it's returning the line that matches plus the two after it. The `-B` flag returns that number of lines before the match, so that's returning the fastq header.

---

**EXERCISE**

Search for the sequence 'GATTTTTACA' (GATTACA with 5 T's instead of 2) in ctl1.fastq file and in the output have the sequence name, e.g.:

```
@SRR1145047.5880759
AAGCTAAAAAAAAAATGGATGTTTCAGTTAAATGTTTTAAAGAGGTACAGATTTTTACAAGGACATAATATAAG
```

Next, search for that sequence in all the FASTQ files.

---

### Redirection & Pipes

We're excited we have all these sequences that we care about that we just got from the FASTQ files. Perhaps that was some really important motif that is going to help us answer some important question. But all those sequences just went whizzing by with grep. How can we capture them?

We can do that with something called "redirection". The idea is that we're redirecting the output to the terminal (all the stuff that went whizzing by) to something else. In this case, we want to write it to a file, so that we can look at it later.

We do this with: `>`.

Let's try it out and put all the sequences that contain 'GATTACA' from all the files in to another file called "gattaca-reads.txt".

```
grep GATTACA *.fastq > gattaca-reads.txt
```

The prompt should sit there a little bit, and then it should look like nothing happened. But type `ls`. You should have a new file called "gattaca-reads.txt". Take a look at it and see if it has what you think it should.

```
ls
less gattaca-reads.txt
```

There's another useful redirection command that we're going to show, and that's called the pipe: `|`. It's probably not a key on your keyboard you use very much. What `|` does is take the output that scrolling by on the terminal and then can run it through another command. When it was all whizzing by before, we wished we could just slow it down and look at it, like we can with `less`. Well it turns out that we can! We pipe the `grep` command through `less`

The pipe '|' takes the output of the first thing and then puts it in to the second part

```
grep GATTACA *.fastq | less
```

Now we can use the arrows to scroll up and down and use `q` to get out.

There's another command called `wc`. `wc` stands for`word count`. It counts the number of lines or characters. So, we can use it to count the number of lines we're getting back from our `grep` command. And that will tell us how many sequences we're finding with that motif across all the files.

```
grep GATTACA *.fastq | wc
```

That tells us the number of lines, words and characters in the file. If we just want the number of lines, we can use the `-l` flag for `lines`.

```
grep GATTACA *.fastq | wc -l
```

You could have piped the output to a file like you did the first time and ran `wc` on that file:

```
wc gattaca-reads.txt
```

But using the pipe you didn't need to create the intermediate file. Pipes are really powerful for stringing together these different commands, so you can do whatever you need to do.

The philosophy behind these command line programs is that none of them really do anything all that impressive. BUT when you start chaining them together, you can do some really powerful things really efficiently. If you want to be proficient at using the shell, you must learn to become proficient with the pipe and redirection operators.

For example, draw a cow saying the largest word in the english language that contains the word "purple":

```bash
# Read in the system's dictionary.
cat /usr/share/dict/words
# grep for "purple"
cat /usr/share/dict/words | grep purple
# count letters in each word
cat /usr/share/dict/words | grep purple | awk '{print length($1), $1}'
# numeric sort
cat /usr/share/dict/words | grep purple | awk '{print length($1), $1}' | sort -n
# take the last line
cat /usr/share/dict/words | grep purple | awk '{print length($1), $1}' | sort -n | tail -n 1
# cut out the second field (-f 2) if the file were delimited by a space (-d " ")
cat /usr/share/dict/words | grep purple | awk '{print length($1), $1}' | sort -n | tail -n 1 | cut -d " " -f 2
# make the cow say it
cat /usr/share/dict/words | grep purple | awk '{print length($1), $1}' | sort -n | tail -n 1 | cut -d " " -f 2 | cowsay
```

This program piped together 7 UNIX programs (`cat`, `grep`, `awk`, `sort`, `tail`, `cut`, `cowsay`), each of which does a very small job on its own, and strung them together to do something pretty powerful (as silly as it may be).

```
 _____________
< unimpurpled >
 -------------
        \   ^__^
         \  (oo)\_______
            (__)\       )\/\
                ||----w |
                ||     ||
```

### Creating, moving, copying, and removing

Now we can move around in the file structure, look at files, search files, redirect. But what if we want to do normal things like copy files or move them around or get rid of them. Sure, if we had a desktop, we could do most of these things without the command line, but what fun would that be?! Besides it's usually faster to do it at the command line, or you'll be on a remote server like Amazon where you won't have another option.

The `coldata.csv` file tells us which sample names are which treatment (in this example the filenames are pretty informative, but when they come off the sequencer usually they won't be). This is a really important file, so we want to make a copy so we don't lose it.

Lets copy the file using the `cp` command. The `cp` command backs up the file. Navigate to the `rnaseq-1day/data` directory and enter:

```
cp coldata.csv coldata.csv-backup
```

Now `coldata.csv-backup` has been created as a copy of `coldata.csv`.

Let's make a `backup` directory where we can put this file.

The `mkdir` command is used to make a directory. Just enter `mkdir` followed by a space, then the directory name.

```
mkdir backup
```

We can now move our backed up file in to this directory. We can move files around using the command `mv`. Enter this command:

```
mv coldata.csv-backup backup
```

This moves `coldata.csv-backup` into the directory `backup/` or the full path would be `~/workshops/lessons/rnaseq-1day/data/backup/coldata.csv-backup`

The `mv` command is also how you rename files. Since this file is important, let's rename it:

```
ls
mv coldata.csv coldata-IMPORTANT.csv
ls
```

Now let's delete the backup copy:

```
ls
ls backup
rm backup/coldata.csv-backup
ls
ls backup
```

The `rm` file removes the file. Be careful with this command. It doesn't just nicely put the files in the Trash. They're really gone.

---

**EXERCISE**

Do the following:

1.	Rename the `coldata-IMPORTANT.csv` file back to `coldata.csv`.
2.	Create a new directory in the rnaseq-1day directory called `new`.
3.	Then, copy the `coldata.csv` file into `new`

---

By default, `rm`, will NOT delete directories. You can tell `rm` to delete a directory and everything in it *recursively* using the `-r` option. Let's delete that `new` directory we just made. Enter the following command:

```
rm -r new
```

Be careful with this.

## Writing files

We've been able to do a lot of work with files that already exist, but what if we want to write our own files. Obviously, we're not going to type in a FASTA file, but there are a lot of reasons we'll want to edit a file.

To write in files, we're going to use the program `nano`. We're going to create a file that contains the favorite grep command so you can remember it for later. We'll name this file 'awesome.sh'.

```
nano awesome.sh
```

Enter the following into the file:

```
grep -B 1 -A 2 GATTACA *.fastq
```

At the bottom of nano, you see the "^X Exit". That means that we use Ctrl-X to exit. Type `Ctrl-X`. It will ask if you want to save it. Type `y` for yes. Then it asks if you want that file name. Hit 'Enter'.


Now you've written a file. You can take a look at it with less or cat, or open it up again and edit it.

---

**EXERCISE**

Open 'awesome.sh' and add "echo AWESOME!" (no quotes) after the grep command and save the file.

We're going to come back and use this file in just a bit.

---

### Running programs

Commands like `ls`, `rm`, `echo`, and `cd` are just ordinary programs on the computer. A program is just a file that you can *execute*. The program `which` tells you the location of a particular program. For example:

```
which ls
```

Will return "/bin/ls". Thus, we can see that `ls` is a program that sits inside of the `/bin` directory. Now enter:

```
which find
```

You will see that `find` is a program that sits inside of the`/usr/bin` directory.

So ... when we enter a program name, like `ls`, and hit enter, how does the shell know where to look for that program? How does it know to run `/bin/ls` when we enter `ls`. The answer is that when we enter a program name and hit enter, there are a few standard places that the shell automatically looks. If it can't find the program in any of those places, it will print an error saying "command not found". Enter the command:

```
echo $PATH
```

This will print out the value of the `PATH` environment variable. Notice that a list of directories, separated by colon characters, is listed. These are the places the shell looks for programs to run. If your program is not in this list, then an error is printed. The shell ONLY checks in the places listed in the `PATH` environment variable.

Remember that file where we wrote our favorite grep command in there? Since we like it so much, we might want to run it again, or even all the time. Instead of writing it out every time, we can just run it as a script. Let's try to run that script:

```
awesome.sh
```

You should get an error saying that hello.sh cannot be found. That is because the directory `~/workshops/lessons/rnaseq-1day/data` is not in the`PATH`. You can try again to run the `awesome.sh` program by entering:

```
./awesome.sh
```

Alas, we get `-bash: ./awesome.sh: Permission denied`. This is because we haven't told the computer that it's a program. To do that we have to make it 'executable'. We do this by changing its mode. The command for that is `chmod` - change mode. We're going to change the mode of this file, so that it's executable and the computer knows it's OK to run it as a program.

To run a program, you have to set the right permissions, make it executable rather than just a text file.

```
chmod +x awesome.sh
```

Now we can run the program

```
./awesome.sh
```

Now you should have seen some output, and of course, it's AWESOME! Congratulations, you just created your first shell script! You're set to rule the world.


## Where can I learn more about the shell?

-	Software Carpentry tutorial - [The Unix shell](http://software-carpentry.org/v4/shell-genomics/index.html)
-	The shell handout - [Command Reference](http://files.fosswire.com/2007/08/fwunixref.pdf)
-	[explainshell.com](http://explainshell.com)
-	<http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html>
- <http://stackoverflow.com/>
-	Google - if you don't know how to do something, try Googling it. Other people have probably had the same question.
-	Learn by doing. There's no real other way to learn this than by trying it out. Write your next paper in nano, open pdfs from the command line, automate something you don't really need to automate or actually do.
