---
layout: page
---

# RNA-seq data analysis: reads to counts

In this section of the workshop we will start with raw RNA-seq reads and do everything that we need to do to get a //count matrix// - a table or "spreadsheet" containing the the number of reads mapping to each gene for each sample.

## Objectives

- Quality assessment for raw sequencing read
- Quality control to trim low-quality sequence data from the ends of reads
- Align reads to the genome for all samples
- Count the number of reads mapping to each gene for each sample

## Data management

First, let's extract the fastq files if we haven't done so already:

```
cd workshops/lessons/rnaseq-1day/data
ll
gunzip *.fastq.gz
ll
```

Next, let's run some quality control analysis. Google search "FastQC" to find the [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) package from Simon Andrews. It's a program that pretty much everyone uses for QC analysis that can be run at the command line on a batch of FASTQ files.

First let's get some help on FastQC. Most bioinformatics programs have a `-h` or `--help` option that gives you some very basic documentation on how to run it. These aren't the same as built-in UNIX man pages. For more in-depth documentation you need to read the paper or visit the tool's website to download a manual. For instance, take a look at the [online FastQC documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) and look specifically at the differences between [good Illumina data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and [bad Illumina data](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).

```
fastqc -h
fastqc *.fastq --threads 4 --outdir .
```

Open up your SFTP client and enter your IP address and credentials. Transfer this data over to your computer to view it. You'll see how the last 5 bases of each read is pretty terrible quality, so let's trim those reads. Google "FASTX Toolkit" to find [Greg Hannon's FASTX-Toolkit software](http://hannonlab.cshl.edu/fastx_toolkit/). This lets you do lots of things for manipulating FASTQ and FASTA files. We're just going to use a feature to trim the end of the reads. You can see which tool you need by looking at the [online documentation](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html).

First, let's run `fastx_trimmer` to get some help.

```
fastx_trimmer -h
```

It tells us we can use the `-t` option followed by a number to trim that number of nucleotides off the end of the read. It tells us that the default way to use the program is to feed data off the `STDIN`, which is UNIX-speak for piping in data. It's default output is the `STDOUT`, which is UNIX-speak for "printing" the data to the screen. Let's see what that means:

```
head -n 8 ctl1.fastq
head -n 8 ctl1.fastq | fastx_trimmer -t 5
```

We get some error about an invalid quality score value. Let's learn a little bit about quality scores. Google FASTQ format and click the Wikipedia article, and scroll down to the [section on quality](http://en.wikipedia.org/wiki/FASTQ_format#Quality). Next scroll down to encoding. Quality scores are encoded in a FASTQ file different ways. You just have to know that this is quality score encoding is Sanger, or Phred+33. There are other tools out there that can tell you what kind of quality score you have, and convert between different quality score encodings. But Sanger/Phred+33 is the most common. Unfortunately `fastx_trimmer` doesn't seem to recognize this encoding. Let's Google that error message specifically. Google "fastx_trimmer invalid quality score value" and you'll arrive at a [SEQAnswers](http://seqanswers.com/forums/showthread.php?t=7596) thread that tells you the answer. (Look who provided that answer, by the way, Simon Andrews, the author of FASTQC!). You can also go to the [FASTX-Toolkit homepage](http://hannonlab.cshl.edu/fastx_toolkit/index.html), and if you scroll down a bit, you'll see that back in 2009 this was added as an option, but it's still not documented in the command line help.

I show you this because this is what doing bioinformatics is like in the real world. You think you have the right tool for the job, you try to run it like the documentation says to run it, and it doesn't work. You Google around to try to find the solution to the problem or better documentation. SEQAnswers is a great resource, and so is Biostars.

Let's try that again:

```
head -n 8 ctl1.fastq
head -n 8 ctl1.fastq | fastx_trimmer -t 5 -Q33
```

Look at that! `fastx_trimmer` trimmed 5 bases off the end of two reads coming in from a pipe and spit it back out to the screen. We could do this to the whole file if we did something like `cat ctl1.fastq | fastx_trimmer ... > trimmed_ctl1.fastq`. Or we could use the command line options to take input files and write output files.

```
fastx_trimmer -h
fastx_trimmer -Q33 -t 5 -i ctl1.fastq -o trimmed_ctl1.fastq
```

Let's take a look at those files to make sure it worked right. Also, let's do a word count as a sanity check to make sure we didn't lose any reads or anything.

```
head -n 4 ctl1.fastq trimmed_ctl1.fastq
wc ctl1.fastq trimmed_ctl1.fastq
```

---

**EXERCISE**

Run `ls` on the current directory. We need to clean things up before they get out of hand.

- Remove the trimmed file you just created
- Make a directory called "QC"
- Move all the fastqc output directories and zip files into the QC directory (hint, you can do this all with a wildcard and a single `mv` command)
- Use `find ... | parallel ...` to run `wc` on all the fastq files in parallel

---




```bash
# extract

# fastqc
fastqc -h
fastqc *.fastq --outdir .

# trim
fastx_trimmer -h
# Q33 undocumented, try heading the file
mkdir Untrimmed
mv *.fastq Untrimmed
cd Untrimmed
find *.fastq | parallel --dry-run fastx_trimmer -t 5 -Q33 -i {} -o ../{}
cd ..
head ctl1.fastq
head Untrimmed/ctl1.fastq

# get chromosome 4 fasta
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr4.fa.gz
# google "ensembl download fasta"
mkdir Index
cd Index
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.4.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.chromosome.4.fa.gz
mv Homo_sapiens.GRCh37.75.dna.chromosome.4.fa chr4.fa

# create index
bowtie-build chr4.fa chr4

# While that's building, let's make annotation
cd Annotation
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.4.fa.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

# create index
bowtie-build chr4.fa chr4

# map
find *.fastq | parallel --dry-run tophat --bowtie1 --no-coverage-search -o {}_tophat Index/chr4 {}
more */align_summary.txt

# count
featureCounts -a Annotation/Homo_sapiens.GRCh37.75.gtf -o counts.txt -t exon -g gene_name */accepted_hits.bam
```
