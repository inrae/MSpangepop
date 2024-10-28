"""
Author: Lucien Piat
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak

Description:
This script parses the header of a VCF file, extracts information sections,
and sorts contigs (genomic references) based on priority:
  - Contigs starting with "." have the highest priority.
  - Numeric contigs (e.g., chromosomes) are sorted in ascending order.
  - Alphabetic contigs are sorted lexicographically.
The script saves the organized header to an output file.

This is done beacause VG wants super well organised vcf files

Usage:
python script_name.py input_file output_file
"""

import argparse
import os

def parse_vcf_header(vcf_header):
    """
    Parses the VCF header to extract general info, contigs, and sorts contigs.

    :param vcf_header: str, the VCF header as a string.
    :return: list, containing combined lines of parsed results.
    """

    lines = vcf_header.strip().split('\n')

    first_part = []
    contigs = []
    second_part = []
    
    start = True
    for line in lines:
        if line.startswith("##contig"):
            contigs.append(line)
            start = False 
        elif start:
            first_part.append(line)
        else:
            second_part.append(line)
    
    def sort_key(contig_line):
        """
        Sorting key function to prioritize and order contigs:
        - Highest priority for contigs starting with ".".
        - Numerical sorting for digit-only IDs (e.g., chromosomes).
        - Lexicographic sorting for alphabetic IDs.
        """
        try:
            contig_id = contig_line.split('<ID=')[1].split(',')[0]
        except IndexError:
            contig_id = ""
        
        if contig_id.startswith("."):
            return (0, contig_id)
        elif contig_id.isdigit():
            return (2, int(contig_id))
        else:
            return (1, contig_id)

    sorted_contigs = sorted(contigs, key=sort_key)
    
    combined_lines = first_part + sorted_contigs + second_part
    return combined_lines

def main(input_file, output_file):
    """
    Main function to process VCF header from input and save sorted output.

    :param input_file: str, path to input file with VCF header.
    :param output_file: str, path to output file for sorted header.
    """

    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        return
    if os.path.getsize(input_file) == 0:
        print(f"Error: Input file '{input_file}' is empty.")
        return

    with open(input_file, 'r') as infile:
        vcf_header = infile.read()
        
    parsed_result = parse_vcf_header(vcf_header)

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(parsed_result) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse VCF header and sort contigs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", help="Path to the input VCF header file.")
    parser.add_argument("output_file", help="Path to the output file for results.")

    args = parser.parse_args()
    
    main(args.input_file, args.output_file)
