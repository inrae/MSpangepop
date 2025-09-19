"""
Author: Lucien Piat  
Institution: INRAe  
Project: PangenOak

Description:  
This script checks whether a given FASTA file is properly formatted and meets the minimum number of required contigs.  
If validation is successful, it generates a .fai index file (in the samtools faidx format) at the specified output path.  
If the FASTA file is malformed or does not meet the minimum requirements, the script will terminate with an error.

Arguments:
  --fasta_file     Path to the input FASTA file.
  --output_file    Path where the output .fai file will be written.
  --min_contigs    Minimum required number of contigs (default: 0).

Note:
  This script is designed to prevent the workflow from running simulations on invalid or incomplete FASTA files.
"""

import argparse
from pathlib import Path
from Bio import SeqIO
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning

def validate_fasta(fasta_records, min_contigs):
    """
    Validates FASTA records for format and uniqueness of contig names.

    Parameters:
        fasta_records (list): List of SeqRecord objects.
        min_contigs (int): Minimum required contig count.

    Returns:
        int: The number of contigs.

    Raises:
        ValueError: If contig names are not unique or minimum contig count is not met.
    """
    contig_names = set()
    
    for record in fasta_records:
        if record.id in contig_names:
            raise MSerror(f"Duplicate contig name found: {record.id}, contig names must be unique")
        contig_names.add(record.id)

    contig_count = len(fasta_records)

    if contig_count < min_contigs:
        raise MSerror(f"Error: The FASTA file contains only {contig_count} contigs, but at least {min_contigs} are required.")
    if contig_count > min_contigs: 
        if min_contigs == 1 :
            MSwarning(f"The simulation wil run on the firt contig (increase chr_n to {contig_count} for full simulation)")
        else :
            MSwarning(f"The simulation wil run on the {min_contigs} frist contigs (increase chr_n to {contig_count} for full simulation)")
        

def generate_fai(fasta_records, output_file):
    """
    Generates a simplified .fai index from parsed FASTA records.

    Parameters:
        fasta_records (list of SeqRecord): Parsed FASTA sequences.
        output_file (str): Path to output .fai file.

    Note:
        This version does not include real byte offsets.
        Fields: contig_name, length, dummy_offset, line_bases, line_width
    """
    with open(output_file, 'w') as fai:
        dummy_offset = 0
        default_line_bases = 60
        default_line_width = 61  # Assuming \n line ending
        for record in fasta_records:
            contig_name = record.id
            seq_len = len(record.seq)
            fai.write(f"{contig_name}\t{seq_len}\t{dummy_offset}\t{default_line_bases}\t{default_line_width}\n")
            dummy_offset += seq_len  # This is not real offset but keeps uniqueness


def main():
    parser = argparse.ArgumentParser(description="Check FASTA format and count contigs.")
    parser.add_argument("--fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("--output_file", help="Path to the output FAI file to be generated.")
    parser.add_argument("--min_contigs", type=int, default=0, help="Minimum required number of contigs (default: 0).")

    args = parser.parse_args()

    try:
        # Read and validate FASTA
        fasta_records = MSpangepopDataHandler.read_fasta(args.fasta_file)
        validate_fasta(fasta_records, args.min_contigs)

        # Generate .fai file at the location of output_file
        generate_fai(fasta_records, args.output_file)


        MSsuccess(f"FASTA file is valid. FAI index created")

    except MSerror as e:
        raise e


if __name__ == "__main__":
    main()