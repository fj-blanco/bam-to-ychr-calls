#!/bin/bash

read -p "Do you want to install SAMtools and BCFtools using apt? (y/n): " install_via_apt

sudo apt-get update
sudo apt-get upgrade

if [ "$install_via_apt" == "y" ]; then
    sudo apt-get install -y samtools
    sudo apt-get install -y bcftools
else
    sudo apt-get install gcc make libncurses5-dev liblzma-dev zlib1g-dev libbz2-dev libcurl4-openssl-dev
    # Install SAMtools from source
    wget -O /tmp/samtools-1.16.1.tar.bz2 https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
    tar -vxjf /tmp/samtools-1.16.1.tar.bz2 -C /tmp
    cd /tmp/samtools-1.16.1
    ./configure --prefix=/usr/local
    make
    sudo make install

    # Install BCFtools from source
    wget -O /tmp/bcftools-1.16.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2
    tar -vxjf /tmp/bcftools-1.16.tar.bz2 -C /tmp
    cd /tmp/bcftools-1.16
    ./configure --prefix=/usr/local
    make
    sudo make install
fi

cd bam-to-ychr-calls
mkdir -p ./genomes/ ./SNP_annotations_hg38/ ./references/

# Download SNP annotations if they don't exist
if [ ! -f "./SNP_annotations_hg38/yleaf_new_positions.txt" ]; then
    wget -O ./SNP_annotations_hg38/yleaf_new_positions.txt https://github.com/genid/Yleaf/raw/master/yleaf/data/hg38/new_positions.txt
fi
if [ ! -f "./SNP_annotations_hg38/ybrowse_snps_hg38.csv" ]; then
    wget -O ./SNP_annotations_hg38/ybrowse_snps_hg38.csv http://ybrowse.org/gbrowse2/gff/snps_hg38.csv
fi

# Download the humanG1Kv37 reference if it doesn't exist and unzip it
if [ ! -f "./references/human_g1k_v37.fasta.gz" ]; then
    wget -O ./references/human_g1k_v37.fasta.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
    gzip -d ./references/human_g1k_v37.fasta.gz
fi