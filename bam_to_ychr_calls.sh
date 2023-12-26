#!/bin/bash

# Default values
output_dir="./results/"
all_calls="false"
quality_threshold=""

# Mandatory flags
sample=""
reference=""

# Parse command-line arguments
while getopts "s:r:o:aq:" flag
do
    case "${flag}" in
        a) all_calls="true";;
        s) sample=${OPTARG};;
        r) reference=${OPTARG};;
        o) output_dir=${OPTARG};;
        q) quality_threshold=${OPTARG};;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$sample" ] || [ -z "$reference" ]; then
    echo "Error: Sample and reference are mandatory. Usage:"
    echo "./script.sh -s sample_name -r reference_file_path [-o output_directory_path]"
    exit 1
fi

# Extract the sample name from the sample path
sample_path=$(realpath "${sample}")
sample_name=$(basename "${sample_path}" .bam)

echo "Sample name: ${sample_name}"

# Ensure that the output directory ends with a slash
[[ "${output_dir}" != */ ]] && output_dir="${output_dir}/"
mkdir -p ${output_dir}
echo "Saving results to ${output_dir}"

# Paths
bam_dir="${sample}"
bam_dir_reidx="${output_dir}${sample_name}_reidx.bam"
bam_Y_dir="${output_dir}${sample_name}_Y.bam"
bcf_dir="${output_dir}${sample_name}_Y.bcf"
vcf_dir="${output_dir}${sample_name}_Y.vcf"

# Function to filter and rename chromosomes
filter_chrY() {
    samtools view -b ${bam_dir} Y > ${bam_Y_dir} 2> /dev/null
    if [ $(samtools view -c ${bam_Y_dir}) -gt 0 ]; then
        echo "Y successfully filtered from .bam file."
    else
        samtools view -b ${bam_dir} chrY > ${bam_Y_dir} 2> /dev/null
        if [ $(samtools view -c ${bam_Y_dir}) -gt 0 ]; then
            echo "chrY successfully filtered from .bam file."
            # Rename chromosomes to remove "chr" prefix
            samtools view -H ${bam_Y_dir} | sed -e 's/SN:chr/SN:/g' | samtools reheader - ${bam_Y_dir} > ${output_dir}${sample_name}_Y_renamed.bam
            bam_Y_dir="${output_dir}${sample_name}_Y_renamed.bam"

        else
            echo "Failed to filter Y chromosome."
            exit 1
        fi
    fi
}

index_bam() {
    # Sort the .bam file
    samtools sort ${bam_dir} -o ${bam_dir_reidx}
    # Index the .bam file
    samtools index ${bam_dir_reidx}
}

check_index() {
    if [ ! -f "${bam_dir}.bai" ]; then
        echo "Index file not found for ${bam_dir}."
        read -p "Do you want to index ${bam_dir} using 'samtools index' (it may take some minutes...)? (y/n): " index_choice
        if [ "$index_choice" == "y" ]; then
            index_bam
        else
            echo "Please ensure that the index ${bam_dir}.bai exists and run the program again"
            exit 1
        fi
    fi
}

# Check index for the original BAM file
check_index

# Filter chrY
filter_chrY

# Variant calling (retaining positions with the ancestral genotype)
if [ "$all_calls" == "true" ]; then
    echo "Variant calling, all calls..."
    bcftools mpileup -Ou -f ${reference} ${bam_Y_dir} | bcftools call -m -Ob -o ${bcf_dir}
else
    echo "Variant calling, variants only..."
    bcftools mpileup -Ou -f ${reference} ${bam_Y_dir} | bcftools call -mv -Ob -o ${bcf_dir}
fi

# Converting the .bcf file to .vcf
bcftools view ${bcf_dir} > ${vcf_dir}

# Extracting the .bam header for exploration (if necessary)
#samtools view -H ${bam_Y_dir} > ${output_dir}${sample_name}_bam_header.txt
echo "Variant calling done."

echo "Running vcf annotation..."
# Building Python command
python_cmd="python ./vcf_processing.py -v ${vcf_dir}"

# Append quality threshold to the command if provided
if [ -n "$quality_threshold" ]; then
    python_cmd+=" -q ${quality_threshold}"
fi

# Append all_calls flag to the command if true
if [ "$all_calls" == "true" ]; then
    python_cmd+=" -a"
fi

# Execute Python command
$python_cmd
echo "Done."