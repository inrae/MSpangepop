"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Simple script to check is the fasta file is corectely formated
--fasta_file Path to the input FASTA file.
--output_file Path to the output file where contig count will be stored.
--min_contigs Minimum required number of contigs (default: 0).

This script will stop the workflow if a fasta file is missing chromosomes or sequences
This will avoid runing the simmulation for nothing
"""

import argparse
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
        

def main():
    parser = argparse.ArgumentParser(description="Check FASTA format and count contigs.")
    parser.add_argument("--fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("--output_file", help="Path to the output file where contig count will be stored.")
    parser.add_argument("--min_contigs", type=int, default=0, help="Minimum required number of contigs (default: 0).")

    args = parser.parse_args()
    MScompute("Validating given fasta file")
    try:
        fasta_records = MSpangepopDataHandler.read_fasta(args.fasta_file)
        validate_fasta(fasta_records, args.min_contigs)

        with open(args.output_file, 'w') as out:
            out.write("")

        MSsuccess("FASTA file is correctly formatted")

    except MSerror as e:
        raise e

if __name__ == "__main__":
    main()