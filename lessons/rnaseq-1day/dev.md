# Notes on developing RNA-seq workshop

Prelim agenda:

0800-1000: Introduction to UNIX shell  
1015-1200: QC, alignment, and expression quantitation  
1300-1500: Introduction to R  
1515-1700: QC and differential expression with R


Intro - will need "remedial" linux/R session at beginning?

First download some data, analyze it, find interesting genes/regions, extract very small fastq from bam files so can do alignment in reasonable time frame. Also, get data only from a single chromosome, so can index that chromosome for aligner so can do this with small amount of RAM on a VM running on an average laptop.

To do:
* convert analysis .R script to .Rmd
* write readme
* adapt intro linux/R materials from SWC or other materials
* load everything on to VM and test
* figure out registration w/ eventbrite
* pick a date


## VM

Stuff need install on VM

Before installing guest additions

```bash
sudo apt-get install gcc make ruby curl git vim cowsay wamerican wamerican-huge wamerican-large
# sudo apt-get install dkms build-essential linux-headers-generic # linux-headers-`uname -r`
# sudo apt-get install virtualbox-guest-x11
```

For my own sanity:

```bash
echo 'alias l="ls -lhGgop --color=always"' >> ~/.bashrc
```

https://github.com/Homebrew/homebrew-science

* bowtie
* tophat
* STAR
* R
* RStudio
* fastqc
* fastx toolkit



Download human genome chromosome 4 from ucsc in case poor internet connection
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr4.fa.gz


etc.

## Orig data

Paper: http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0093338

GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54413

Put these files into a text file called `sra-files.txt`

```
# UVB
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450265/SRR1145042/SRR1145042.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450266/SRR1145043/SRR1145043.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450267/SRR1145044/SRR1145044.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450268/SRR1145045/SRR1145045.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450269/SRR1145046/SRR1145046.sra

# Control
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450270/SRR1145047/SRR1145047.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450271/SRR1145048/SRR1145048.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450272/SRR1145049/SRR1145049.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450273/SRR1145050/SRR1145050.sra
```

Then use wget to get all the files

```bash
wget sra-files.txt
```

Now use the sra-toolkit to unpack all the files.


```bash
find *sra | parallel fastq-dump {}
```

Use a perl script to convert all to Illumina 1.9

```bash
find *fastq | parallel "illumina15-to-sanger.pl {} > {}.fq"
```

FASTQC on all to make sure the quality encoding worked well.

```bash
find *.f*q | parallel fastqc {} --outdir .
```

```bash
mv SRR1145042.fastq.fq uvb1.fq
mv SRR1145043.fastq.fq uvb2.fq
mv SRR1145044.fastq.fq uvb3.fq
mv SRR1145045.fastq.fq uvb4.fq
mv SRR1145046.fastq.fq uvb5.fq

mv SRR1145047.fastq.fq ctl1.fq
mv SRR1145048.fastq.fq ctl2.fq
mv SRR1145049.fastq.fq ctl3.fq
mv SRR1145050.fastq.fq ctl4.fq


## Alignment
NCPUS=12
MEM="40GB"
GENOMEDIR="/home/sdt5z/genomes/star/hg19/"
QSUB="qsub -l select=1:ncpus=$NCPUS:mem=$MEM,walltime=24:00:00 -q uvabx -W group_list=uvabx -V -j oe -m bae -M vustephen+fir@gmail.com"
find `pwd` -name "*.fq" | sed 's/\.fq$//g' | sort | xargs -i echo $QSUB -- `which time` `which STAR` --genomeDir $GENOMEDIR --runThreadN $NCPUS --outFileNamePrefix {}. --readFilesIn {}.fq > runstar.sh

## Count
# featureCounts -a ~/genomes/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o counts.txt -T 12 -t exon -g gene_id *.sam

## Extract chr4:60000000-160000000
echo -e "chr4\t60000000\t160000000" > roi.bed
find *.sam | sed 's/.Aligned.out.sam//g' | sort | parallel --dry-run 'samtools view -Sb {}.Aligned.out.sam | bedtools intersect -abam - -b roi.bed | bedtools bamtofastq -i - -fq {}.fastq'
```


## New small data

```bash
# extract
gunzip *.fastq.gz

# fastqc
find *.fastq | parallel --dry-run fastqc {} --outdir .

# trim
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
