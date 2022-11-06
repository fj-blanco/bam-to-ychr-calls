
## Description

This code performs variant calling with [BCFtools](https://samtools.github.io/bcftools/) from a .bam genome in hg19 coordinates, keeping the ancestral positions. Then it annotates the vcf using the YLeaf and Ybrowse lists of annotated SNPs in hg38 coordinates. The change of coordinate is done with [pyliftover](https://pypi.org/project/pyliftover/). The result is a .csv with annotated SNPs sorted by Phred quality score.

In the process we annotate the regions not suitable for phylogeny [as per YSEQ](https://www.yseq.net/product_info.php?products_id=108&osCsid=a46df681c44538157cf8939e4aeef532):

- **Pseudo Autosomal Region 1** (PAR1): chrY:1..2781479
- **Synthetic Assembled Centromeric Region** (CEN): chrY:10072350..11686750
- **DYZ19 125 bp repeat region**: chrY:20054914..20351054 (DYZ19 125 bp repeat region)
- **Post Palindromic Region**: chrY:26637971..26673210 
- **Pseudo Autosomal Region 2**, (PAR2): chrY:56887903..57217415

It also annotates mutations prone to aDNA damage ["C to T", "A to T", "G to C", "A to G"]

## Requirements

SAMtools  
BCFtools

## Setup

The following installation has been tested in Debian GNU/Linux 11.

## Installing SAMtools and BCFtools:
You may want to check this useful links: [link1](https://www.biostars.org/p/328831/) and [link2](https://www.biostars.org/p/328831/).

You may want to check if there are [newer SAMtools & BCFtools releases](https://www.htslib.org/download/). To install another version just make the relevant replacements in the code below.

### Update, upgrade and install dependencies:
```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install gcc
sudo apt-get install make
sudo apt-get install libncurses5-dev
sudo apt-get install liblzma-dev
sudo apt-get install libbz2-dev
apt-get install libcurl4-openssl-dev
```
### Install SAMtools:

```
cd /usr/bin
sudo wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
sudo tar -vxjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
sudo make
sudo make install
export PATH="$PATH:/usr/bin/samtools-1.16.1"
source ~/.profile
```
### Install BCFtools:

```
cd ..
sudo wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
sudo tar -vxjf bcftools-1.16.tar.bz2
cd bcftools-1.16
sudo make
sudo make install
export PATH="$PATH:/usr/bin/bcftools-1.16"
source ~/.profile
```

## Installing the repository
### Clone this repository:
```
git clone https://github.com/fj-blanco/bam-to-ychr-calls
```

### Move to the repository directory and create the necessary folders:
```
cd bam-to-ychr-calls
sudo mkdir ./genomes/
sudo mkdir ./SNP_annotations_hg_38/
sudo mkdir ./references/
```

### Download the SNP annotations:
```
wget -O ./SNP_annotations_hg38/yleaf_new_positions.txt https://github.com/genid/Yleaf/raw/master/yleaf/data/hg38/new_positions.txt
wget -O ./SNP_annotations_hg38/ybrowse_snps_hg38.csv http://ybrowse.org/gbrowse2/gff/snps_hg38.csv
```
### Download the humanG1Kv37 reference and unzip it:
```
wget -O ./references/human_g1k_v37.fasta.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gzip -d ./references/human_g1k_v37.fasta.gz
```

### Download the genome:
As an example we will use cay0007 from the recent paper [A genomic snapshot of demographic and cultural dynamism in Upper Mesopotamia during the Neolithic Transition](https://www.science.org/doi/10.1126/sciadv.abo3609).
```
wget -O ./genomes/sample.bam ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR957/ERR9576228/cay007_merged_220118.hs37d5.cons.90perc.bam
```

### Create a Python virtual environment (optional) and install the dependencies:
```console
python3 -m venv bam-to-ychr-calls-venv
source bam-to-ychr-calls-venv/bin/activate
pip install --no-cache-dir -r requirements.txt
```

### Giving execute permissions and executing the bash script:
```
chmod +x bam_to_vcf.sh
./bam_to_vcf.sh
```
### Run the Python script:
```
python3 ./vcf_processing.py
```