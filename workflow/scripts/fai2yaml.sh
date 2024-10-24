#!/bin/bash

# Author: Lucien Piat
# Date: October 24, 2024
# Project: PangenOak at INRAE

# Check if the input file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fai_file> <output_yaml_file>"
    exit 1
fi

input_fai_file="$1"
output_yaml_file="$2"

# Initialize the output YAML file
echo "chromosomes:" > "$output_yaml_file"

# Read the FAI file and extract chromosome names
while IFS=$'\t' read -r chromosome _; do
    # Append each chromosome name to the YAML file
    echo "  - \"$chromosome\"" >> "$output_yaml_file"
done < "$input_fai_file"

echo "YAML index file generated: $output_yaml_file"
