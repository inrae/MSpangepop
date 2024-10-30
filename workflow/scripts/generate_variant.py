"""
Author: Lucien Piat, based on Sukanya Denni's initial script
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak

This script generates genetic variants based on input data from msprime simulated VCF, FAI, and FASTA files.
It reads variant probabilities and lengths from specified files, creates variants for each 
row in the VCF file, and writes the generated variants to a new VCF file.
"""

import random
import pandas as pd
import argparse
from readfile import read_fasta, read_variant_length_file, read_fai, read_vcf, read_yaml
from variants_class import *

def select_length(df):
    """Select a length from the interval based on cumulative probabilities and return an integer length."""
    rand_val = random.random()
    row = df[df['cumulative_pb'] >= rand_val].iloc[0]
    interval = row['size_interval']
    lower_bound, upper_bound = map(float, interval.strip('[]').split(','))
    return random.randint(int(lower_bound), int(upper_bound))

def retrieve_reference_sequence(chromosome_dict, chrom, pos, length, no_ref=False):
    """Retrieve the reference sequence from a FASTA file."""
    if chrom in chromosome_dict:
        seq_record = chromosome_dict[chrom]
        start = pos - 1 
        if no_ref:
            return str(seq_record.seq[start:start + 1])  # Only need one base for insertions
        end = min(start + length, len(seq_record))
        return str(seq_record.seq[start:end])
    else:
        raise ValueError(f"Chromosome {chrom} not found in FASTA.")

def create_variant(variant_probs, chrom, pos, length_files, chrom_lengths, chromosome_dict, samples):
    """Create a variant with a random length based on probabilities and file data, respecting chromosome length."""
    variant_types = list(variant_probs.keys())
    probabilities = list(variant_probs.values())
    
    selected_variant_type = random.choices(variant_types, weights=probabilities, k=1)[0]
    
    # Dynamically get the class using globals() and selected_variant_type
    variant_class = globals().get(selected_variant_type)
    
    if variant_class is None:
        raise ValueError(f"Variant type '{selected_variant_type}' does not correspond to an implemented class.")
    
    if selected_variant_type == 'SNP':
        variant = variant_class(chrom, pos)
    else:
        length_df = length_files[selected_variant_type]
        length = select_length(length_df)
        
        # Ensure the variant length does not exceed the chromosome length
        max_length = chrom_lengths[chrom] - pos
        if length > max_length:
            length = max_length
        if length == 0:
            length += 1
            
        variant = variant_class(chrom, pos, length)

    # Retrieve the reference sequence and attach it to the variant
    no_ref = (selected_variant_type == 'Insertion')
    variant.reference_seq = retrieve_reference_sequence(chromosome_dict, chrom, pos, variant.length, no_ref)
    
    variant.samples = samples
    variant.compute_alt_seq()
    
    return variant


def main(fai_file, vcf_input_file, vcf_output_file, fasta_file, yaml_file):
    # Load variant probabilities and length files dynamically from YAML input
    variant_probs = read_yaml(yaml_file)

    length_files = {
        'Deletion': read_variant_length_file('simulation_data/size_distribDEL.tsv'),
        'Insertion': read_variant_length_file('simulation_data/size_distribINS.tsv'),
        'Inversion': read_variant_length_file('simulation_data/size_distribINV.tsv'),
        'TandemDuplication': read_variant_length_file('simulation_data/size_distribDUP.tsv'),
        'InvertedTandemDuplication': read_variant_length_file('simulation_data/size_distribDUP.tsv'),
        'Translocation': read_variant_length_file('simulation_data/size_distribINS.tsv'),
        'Transduplication': read_variant_length_file('simulation_data/size_distribINS.tsv'),
        'ReciprocalTranslocation': read_variant_length_file('simulation_data/size_distribINS.tsv')
    }

    chrom_lengths = read_fai(fai_file)
    vcf_data = read_vcf(vcf_input_file)
    chromosome_dict = read_fasta(fasta_file)

    with open(vcf_output_file, 'w') as vcf_output:
        # Generate and write variants for each row in the VCF
        for index, random_row in vcf_data.iterrows():
            chrom = random_row['CHROM']
            pos = random_row['POS']
            
            samples = {col: random_row[col] for col in vcf_data.columns if col.startswith('SAMPLE')}
            
            variant_instance = create_variant(variant_probs, chrom, pos, length_files, chrom_lengths, chromosome_dict, samples)
            vcf_line = variant_instance.vcf_line()

            vcf_output.write(vcf_line + '\n')

if __name__ == "__main__":
    # Set up argument parsing with added YAML file path
    parser = argparse.ArgumentParser(description="Generate variants from VCF and other genomic data.")
    parser.add_argument("--fai", help="Path to the FAI file.")
    parser.add_argument("--vcf", help="Path to the input VCF file.")
    parser.add_argument("--output", help="Path to the output VCF file.")
    parser.add_argument("--fasta", help="Path to the FASTA file.")
    parser.add_argument("--yaml", help="Path to the YAML file with variant probabilities.")

    args = parser.parse_args()
    main(args.fai, args.vcf, args.output, args.fasta, args.yaml)
