#!/bin/bash
# .bam to ychr .vcf:

sample="sample";
samples_directory="./genomes/";
reference="./references/human_g1k_v37.fasta"

bam_dir="${samples_directory}${sample}.bam"
bam_dir_sorted="${samples_directory}${sample}_sorted.bam"
bam_Y_dir="${samples_directory}${sample}_Y.bam"
bcf_dir="${samples_directory}${sample}_Y.bcf"
vcf_dir="${samples_directory}${sample}_Y.vcf"

# Ordering the .bam file:
samtools sort ${bam_dir} -o ${bam_dir_sorted}

# Indexing the .bam file:
samtools index ${bam_dir_sorted}

# We keep just the Y chromosome
samtools view -b ${bam_dir_sorted} Y > ${bam_Y_dir}

# Indexig the Y chromosome .bam file:
samtools index ${bam_Y_dir}

# Variant calling (retaining positions with the ancestral genotype).
# Change to  bcftools call -mv -Ob -o calls.bcf for variants only:
bcftools mpileup -Ou -f ${reference} ${bam_Y_dir} | bcftools call -m -Ob -o ${bcf_dir}

# Converting the .bcf file to .vcf:
bcftools view ${bcf_dir} > ${vcf_dir}

# Extracting the .bam header for exploration (if necessary):
samtools view -H ${bam_Y_dir} > ${samples_directory}${sample}_bam_header.txt