#!/bin/bash

# This script takes a FASTA file as input, splits each contig into individual smaller FASTA files.
#
# Created by Lucien Piat on behalf of INRAe as part of the PangenOak project.
#
# Usage:
#   ./split_fasta.sh <input_fasta> <output_directory> [--help]
#
# Each output file is named after the contig header (without the ">"), with a `.fasta` extension.

# Display help message
show_help() {
    echo "Usage: $0 <input_fasta> <output_directory> [--help]"
    echo
    echo "This script splits a FASTA file into individual contig files."
    echo "Each output file is named after the contig header and saved in the specified output directory."
    echo
    echo "Arguments:"
    echo "  <input_fasta>       Path to the input FASTA file."
    echo "  <output_directory>  Directory to save the individual FASTA files."
    echo
    echo "Options:"
    echo "  --help              Display this help message and exit."
    echo
    exit 0
}

if [[ "$1" == "--help" ]]; then
    show_help
fi

if [[ "$#" -ne 2 ]]; then
    echo "Error: Missing arguments."
    show_help
fi

in="$1"
out_dir="$2"

mkdir -p "$out_dir"

out=""
while read -r line; do
    if [[ "$line" == \>* ]]; then
        # Close previous file descriptor if open
        [ -n "$out" ] && exec 3>&-
        # Prepare the new output file
        out="${line#>}.fasta"
        exec 3> "$out_dir/$out"
        echo "$line" >&3
    else
        # Write sequence data to the current output file
        echo "$line" >&3
    fi
done < "$in"

[ -n "$out" ] && exec 3>&-

echo "FASTA file splitting complete. Files saved in '$out_dir'."