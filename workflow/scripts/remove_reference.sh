#!/bin/bash

# Author: Lucien Piat
# Date: November 24, 2024
# Project: PangenOak at INRAE

# Usage: ./remove_reference.sh input_fasta output_fasta
input_fasta=$1
output_fasta=$2

# Filter the input FASTA file to include only sequences starting with 'tsk'
awk 'BEGIN {RS=">"; ORS=""} /^tsk/ {print ">" $0}' "$input_fasta" > "$output_fasta"
