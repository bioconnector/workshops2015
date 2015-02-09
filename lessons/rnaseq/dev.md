# Notes on developing RNA-seq workshop

## Setting up the machine image

Primarily using AWS image for course. Also creating Lubuntu VirtualBox image as backup to install on local computers just in case we have problems with AWS.

### AWS

First, create a new user **bioinfo** and log in as this user. Also, set the password for the **ubuntu** user.

```bash
# create ubuntu password, set to ubuntu
sudo passwd ubuntu

# create user bioinfo, set password to bioinfo, add to sudoers
sudo adduser bioinfo
sudo adduser bioinfo sudo

# login as bioinfo to do the rest
su - bioinfo
```

Now run the software installation script below. Then come back and lock things down. Once you run this you won't be able to log into the same instance again. Save the image and launch a new instance from the image.

```bash
# Allow password login
sudo vim /etc/ssh/sshd_config # PasswordAuthentication yes
sudo service ssh restart

# Remove ssh host key pairs and authorized keys
sudo shred -u /etc/ssh/*_key /etc/ssh/*_key.pub
sudo find / -name "authorized_keys" -exec rm -f {} \;

# Remove all shell history
history -w
history -c
shred -u ~/.*history
sudo find /root/.*history /home/*/.*history -exec rm -f {} \;
history -w
history -c
```

Now save the machine image.

* **Name:** rnaseq
* **Description:** Ubuntu image for UVA Bioconnector RNA-seq workshop (samtools, bowtie, tophat2, featureCounts, fastx_toolkit, fastqc, etc)


### Virtualbox backup image

Just in case we have a problem with AWS, let's create a VirtualBox image using Lubuntu that looks pretty much the same as the AWS image.

First, install the linux headers and the VirtualBox guest additions. Then install R, RStudio, and select R packages.

```bash
# If using virtualbox and want to install guest additions
sudo apt-get -y install gcc make linux-headers-generic linux-headers-$(uname -r)

# Install guest additions on /media/user/disc
sudo ./VBoxLinuxAdditions.run
sudo reboot

# Install R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
sudo echo "deb http://http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
sudo apt-get -y install r-base r-cran-xml

# Install RStudio
wget http://download1.rstudio.org/rstudio-0.98.1087-amd64.deb
sudo dpkg -i rstudio-0.98.1087-amd64.deb
```

Within R

```r
options("repos" = c(CRAN = "http://cran.rstudio.com/"))
install.packages("calibrate")
install.packages("gplots")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("DESeq2")
```
After finishing setup right before saving the image, clear the history:

```bash
# Remove all shell history
history -w
history -c
shred -u ~/.*history
sudo find /root/.*history /home/*/.*history -exec rm -f {} \;
history -w
history -c
```

### Setup script

Use this on both the AWS or VirtualBox image to:

* Update & upgrade
* Install system software
* Install samtools, fastx-toolkit, and fastqc from apt.
* Download binaries for bowtie2, tophat2, and subread, and add these binaries to path.
* Download and extract the genomedata directory used in the course that contains the chromosome 4 fasta file, GTF file, and bowtie2 indexes.
* Clone the workshop directory.

```bash
############################## install.sh ######################################

# install from apt
sudo apt-get -y update
sudo apt-get -y upgrade
sudo apt-get -y install gcc make ruby curl git vim parallel unzip firefox cowsay wamerican-huge default-jre
sudo apt-get -y install samtools fastx-toolkit fastqc

# Download and extract software manually
cd
mkdir -p bin
cd bin
## bowtie2
wget http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip
unzip bowtie2-2.2.4-linux-x86_64.zip
ln -s bowtie2-2.2.4 bowtie2
## tophat2
wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.13.Linux_x86_64.tar.gz
tar zxvf tophat-2.0.13.Linux_x86_64.tar.gz
ln -s tophat-2.0.13.Linux_x86_64 tophat
## featureCounts
wget http://downloads.sourceforge.net/project/subread/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz
tar zxvf subread-1.4.6-Linux-x86_64.tar.gz
ln -s subread-1.4.6-Linux-x86_64.tar.gz subread

# Change path
cd
echo 'export PATH=$HOME/bin/bowtie2-2.2.4/:$HOME/bin/subread-1.4.6-Linux-x86_64/bin/:$HOME/bin/tophat-2.0.13.Linux_x86_64/:$PATH' >> ~/.bashrc

# Get data here
cd
wget http://people.virginia.edu/~sdt5z/genomedata.tgz .
tar xvf genomedata.tgz
rm genomedata.tgz

# Get workshop data
git clone https://github.com/bioconnector/workshops.git

################################################################################
```

## Work done for setting up

### Create genomedata

```bash
# download and extract genome data:
mkdir -p genomedata
cd genomedata
wget -b ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz
wget -b ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
gunzip *.gz
mv Homo_sapiens.GRCh38.dna.chromosome.4.fa chr4.fa
mv Homo_sapiens.GRCh38.77.gtf genes.gtf
grep ^4 genes.gtf > chr4.gtf
gzip genes.gtf
samtools faidx chr4.fa
bowtie2-build chr4.fa chr4
```

### Original data to small data

This section outlines how I took the original data, did the DGE analysis, then selected just the regions with interesting hits for the example we'll do in class.

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

Rename the files to something meaningful.

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
```

Align, count, run the DGE, then extract regions of interest.

```bash
## Alignment
NCPUS=12
MEM="40GB"
GENOMEDIR="/home/sdt5z/genomes/star/hg19/"
QSUB="qsub -l select=1:ncpus=$NCPUS:mem=$MEM,walltime=24:00:00 -q uvabx -W group_list=uvabx -V -j oe -m bae -M vustephen+fir@gmail.com"
find `pwd` -name "*.fq" | sed 's/\.fq$//g' | sort | xargs -i echo $QSUB -- `which time` `which STAR` --genomeDir $GENOMEDIR --runThreadN $NCPUS --outFileNamePrefix {}. --readFilesIn {}.fq > runstar.sh

## Count
featureCounts -a ~/genomes/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o counts.txt -T 12 -t exon -g gene_id *.sam

## Extract chr4:60000000-160000000
echo -e "chr4\t60000000\t160000000" > roi.bed
find *.sam | sed 's/.Aligned.out.sam//g' | sort | parallel --dry-run 'samtools view -Sb {}.Aligned.out.sam | bedtools intersect -abam - -b roi.bed | bedtools bamtofastq -i - -fq {}.fastq'
```
