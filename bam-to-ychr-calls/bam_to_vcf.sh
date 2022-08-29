#!/bin/bash

# VCF TO 23andMe:

# Ejecutar con (en environment_1):
# ./bam_to_csv.sh
# Ejecutar en el directorio:
# cd ./Documents/variant_calling/ ;

# Modificar la etiqueta de individuo:
sample="I13518";
samples_directory="/home/javi/Documents/genomics/samples/southern_arc/${sample}/";
#reference="hg38_ncbi.fa";
reference="human_g1k_v37.fasta"
#reference="./other_references/hg19_reference.fa";

bam_dir="${samples_directory}${sample}.bam"
bam_dir_sorted="${samples_directory}${sample}_sorted.bam"
bam_Y_dir="${samples_directory}${sample}_Y.bam"
bcf_dir="${samples_directory}${sample}_Y.bcf"
vcf_dir="${samples_directory}${sample}_Y.vcf"

# Ordenamos el bam:
samtools sort ${bam_dir} -o ${bam_dir_sorted}

# Indexamos el BAM:
samtools index ${bam_dir_sorted}

# Filtramos para conservar solo el cromosoma Y:
samtools view -b ${bam_dir_sorted} Y > ${bam_Y_dir}

# Indexamos el bam con el cromosoma Y:
samtools index ${bam_Y_dir}

# Solo variantes: usar el parametro -mv despues de call, como en
# bcftools call -mv -Ob -o calls.bcf
# Para retener posiciones con genotipo ancestral usar el parametro -m despues de call, como en
# bcftools call -m -Ob -o calls.bcf
bcftools mpileup -Ou -f ${reference} ${bam_Y_dir} | bcftools call -m -Ob -o ${bcf_dir}

# Para pasar a vcf:
bcftools view ${bcf_dir} > ${vcf_dir}


# Miramos la cabecera del bam:
#samtools view -H ${bam_Y_dir} > ${samples_directory}${sample}_bam_header.txt