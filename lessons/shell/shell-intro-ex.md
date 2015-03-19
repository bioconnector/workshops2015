# Exercises for shell intro

**Exercise 1**

- Using `cd` and `ls`, go in to the `workshops/shell/data` directory and list its contents.
- How many files, how many directories and how many programs are there?

----

**Exercise 2**

Let's go on a file hunt. Move around in the `shell/data/hidden` directory and try to find the file `youfoundit.txt`.

----

**Exercise 3**

Try finding the `anotherfile.txt` file without changing directories.

----

**Exercise 4**

- List the contents of the `/bin` directory. Do you see anything familiar in there?

----

**Exercise 5**

Do each of the following using a single `ls` command without navigating to a different directory.
  - List all of the files in `/bin` that start with the letter 'c'
  - List all of the files in `/bin` that contain the letter 'a'
  - List all of the files in `/bin` that end with the letter 'o'

----

**Exercise 6**

Find the line number in your history for the last exercise (listing files in `/bin`) and reissue that command.

----

**Exercise 7**

- Print out the contents of the `workshops/lessons/rnaseq/data/coldata.csv` file. What does this file contain?
- Without changing directories, use one short command to print the contents of all of the files in the `~/workshops/posts_/` directory.

----

**Exercise 8**

Search for the sequence "GATTTTTACA" (GATTACA with 5 T's instead of 2) in ctl1.fastq file and in the output have the sequence name. The output should look like this:

```
@SRR1145047.5880759
AAGCTAAAAAAAAAATGGATGTTTCAGTTAAATGTTTTAAAGAGGTACAGATTTTTACAAGGACATAATATAAG
```

Next, search for that sequence in all the FASTQ files.

----

**Exercise 9**

Do the following:

1.  Rename the `coldata-IMPORTANT.csv` file back to `coldata.csv`.
2.  Create a new directory in the rnaseq directory called `new`.
3.  Then, copy the `coldata.csv` file into `new`

----

**Exercise 10**

Open `awesome.sh` and add `echo AWESOME!` after the grep command and save the file.

We're going to come back and use this file in just a bit.

----

**Exercise 11**

1. In the `data` directory, use `nano` to write a script called `quickpeek.sh` that:
    * Runs `head` on all the fastq files in the current directory
    * Runs `wc` on all the fastq files
    * `echo`s "Done!" when finished.
2. Make the program executable.
3. Run the program.

----

**Exercise 12**

1. `cd` into the `data` directory and take a look at what files were created.
2. Open up one of the new files with `less` and use the `/` key to search for the "GATTACA" motif. Does it actually occur on every line?
3. Use `rm` to delete all the files ending with `.gattaca.txt`.
