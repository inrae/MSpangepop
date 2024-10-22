#!/bin/bash

# Author: Lucien Piat
# Date: October 24, 2024
# Project: PangenOak at INRAE

# Input parameters
output_file="$1"  # Output merged VCF file
shift               # Shift to get the remaining arguments as input VCF files

# Check if output file is provided
if [ -z "$output_file" ]; then
    echo "Usage: $0 <output_file> <input_vcf1> <input_vcf2> ..."
    exit 1
fi

# Create the output file
> "$output_file"

# Iterate over all provided VCF files
for vcf_file in "$@"; do
    if [ ! -f "$vcf_file" ]; then
        echo "Warning: File $vcf_file does not exist. Skipping."
        continue
    fi
    # If output file is empty (first file), include the column header
    if [ ! -s "$output_file" ]; then
        grep -E '^#CHROM' "$vcf_file" >> "$output_file"
    fi
    # Append the VCF content without headers (skip lines starting with #)
    grep -v '^#' "$vcf_file" >> "$output_file"
done

echo "Merging completed. Output file: $output_file"
