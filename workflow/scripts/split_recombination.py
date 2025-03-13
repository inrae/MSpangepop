"""
Author: Lucien Piat
Creation: 4 Feb 2025
Updated: 4 Feb 2025
Institution: INRAe
Project: PangenOak
"""

import argparse
import sys
import gzip
from readfile import read_json, read_fasta

def extract_intervals(fasta_sequences, json_data, chromosome_number, output_file):
    """
    Extracts sequences for the specified intervals from the FASTA file and saves them to a bgzipped output file.

    Parameters:
        fasta_sequences (list): List of sequences from the FASTA file.
        json_data (list): List of intervals and mutation data from the JSON file.
        chromosome_number (int): The 1-based index of the chromosome to extract.
        output_file (str): Path to the output file (will be bgzipped).
    """
    chrom_index = chromosome_number - 1  # Convert 1-based index to 0-based
    
    if chrom_index < 0 or chrom_index >= len(fasta_sequences):
        print(f"❌ Error: Chromosome {chromosome_number} not found in FASTA file.")
        sys.exit(1)

    try:
        chrom_seq = fasta_sequences[chrom_index].seq  # Retrieve chromosome sequence
    except AttributeError:
        print(f"❌ Error: Invalid sequence format for chromosome {chromosome_number}.")
        sys.exit(1)

    bgzip_output_file = output_file if output_file.endswith(".gz") else output_file + ".gz"

    try:
        with gzip.open(bgzip_output_file, "wt") as out:  # Open file in text mode for writing
            for entry in json_data:
                try:
                    start, end = map(int, entry["interval"])  # Convert float to int for slicing
                    if start < 0 or end > len(chrom_seq):
                        print(f"⚠️ Warning: Interval {start}-{end} is out of bounds for chromosome {chromosome_number}. Skipping.")
                        continue
                    
                    extracted_seq = chrom_seq[start:end]  # Extract sequence for the interval
                    output_line = f">Chromosome{chromosome_number}_Interval_{start}_{end}\n{extracted_seq}\n"
                    out.write(output_line)  # Write to bgzipped file
                    print(f"✅ MSpangepop -> Saved interval {start}-{end} from chromosome {chromosome_number}")
                except (KeyError, TypeError, ValueError) as e:
                    print(f"⚠️ Warning: Skipping invalid JSON entry {entry}. Error: {e}")
    except IOError as e:
        print(f"❌ Error writing to output file {bgzip_output_file}: {e}")
        sys.exit(1)

def main(json_file, fasta_file, chromosome_number, output_file):
    """
    Reads JSON and FASTA files, extracts sequences for specified intervals, and saves them to a bgzipped file.

    Parameters:
        json_file (str): Path to the JSON file.
        fasta_file (str): Path to the FASTA file.
        chromosome_number (int): Chromosome number (1-based index).
        output_file (str): Path to the output file (will be bgzipped).
    """
    try:
        fasta_sequences = read_fasta(fasta_file)
    except Exception as e:
        print(f"❌ Error reading FASTA file: {e}")
        sys.exit(1)

    try:
        json_data = read_json(json_file)
    except Exception as e:
        print(f"❌ Error reading JSON file: {e}")
        sys.exit(1)

    print(f"🔹 MSpangepop -> Starting to split chr {chromosome_number} between recombination.")
    extract_intervals(fasta_sequences, json_data, chromosome_number, output_file)
    print(f"✅ MSpangepop -> Split successful on chr {chromosome_number}, output saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences for given intervals from a FASTA file and save to a bgzipped output file.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--fasta", required=True, help="Path to the FASTA file.")
    parser.add_argument("--chromosome", type=int, required=True, help="Chromosome number (1-based index).")
    parser.add_argument("--output", required=True, help="Path to output file (will be bgzipped).")

    args = parser.parse_args()
    main(args.json, args.fasta, args.chromosome, args.output)
