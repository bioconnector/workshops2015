# Notes on developing RNA-seq workshop

## Machine image

Ubuntu 14.04 LTS image from AWS

```bash
# install software
sudo apt-get -y update
sudo apt-get -y upgrade
sudo apt-get -y install gcc make ruby curl git vim parallel unzip firefox cowsay wamerican-huge
sudo apt-get -y install samtools fastx-toolkit fastqc

# download and extract manually:
# bowtie2
# tophat2
# featureCounts

# download and extract genome data:
mkdir genomedata
cd genomedata
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip *.gz
mv Homo_sapiens.GRCh38.dna.chromosome.4.fa chr4.fa
mv Homo_sapiens.GRCh38.77.gtf genes.gtf
grep ^r genes.gtf > chr4.gtf
# need to index the fa with samtools faidx and then create bowtie2 indexes.
```

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
