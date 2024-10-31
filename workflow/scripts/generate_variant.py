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

def create_variant(variant_probs, chrom, position, length_files, chrom_lengths, chromosome_dict, samples):
    """Create a variant with a random length based on probabilities and file data, respecting chromosome length."""
    variant_types = list(variant_probs.keys())
    probabilities = list(variant_probs.values())
    
    selected_variant_type = random.choices(variant_types, weights=probabilities, k=1)[0]
    
    # Dynamically get the class
    variant_class = globals().get(selected_variant_type)
    if variant_class is None:
        raise ValueError(f"Variant type '{selected_variant_type}' does not correspond to an implemented class.")
    
    # Instantiate the variant object with relevant parameters
    # TODO, include no_ref as a variant variable
    if selected_variant_type == 'SNP':
        variant = variant_class(chrom, position, chromosome_dict, chrom_lengths)
    else:
        length_df = length_files[selected_variant_type]
        variant = variant_class(chrom, position, chromosome_dict, chrom_lengths, length_df)

    no_ref = (selected_variant_type == 'Insertion')
    variant.reference_seq = variant.retrieve_reference_sequence(no_ref=no_ref)
    
    variant.samples = samples
    variant.compute_alt_seq()
    
    # Special handling for Translocation
    if isinstance(variant, Translocation):
        # Choose a random destination chromosome and position
        dest_chrom = random.choice(list(chrom_lengths.keys()))
        dest_position = random.randint(1, chrom_lengths[dest_chrom])
        variant.set_destination(dest_chrom, dest_position)
        
        variant.compute_alt_seq()

    return variant

def main(fai_file, vcf_input_file, vcf_output_file, fasta_file, yaml_file):
    # Load variant probabilities and length files dynamically from YAML input
    variant_probs = read_yaml(yaml_file)

    # Load length distributions for each variant type from files
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

    # Load chromosome lengths, VCF data, and chromosome dictionary
    chrom_lengths = read_fai(fai_file)
    vcf_data = read_vcf(vcf_input_file)
    chromosome_dict = read_fasta(fasta_file)

    # Open the output VCF file for writing
    with open(vcf_output_file, 'w') as vcf_output:
        
        for row in vcf_data.iterrows():
            chrom = row['CHROM']
            position = row['POS']
            
            # Extract sample data from VCF row
            samples = {col: row[col] for col in vcf_data.columns if col.startswith('SAMPLE')}
            
            # Create a variant instance and generate VCF line output
            variant_instance = create_variant(variant_probs, chrom, position, length_files, chrom_lengths, chromosome_dict, samples)
            vcf_line = variant_instance.vcf_line()
            print(variant_instance.describe())
            vcf_output.write(vcf_line + '\n')

if __name__ == "__main__":
    # Set up argument parsing with YAML file path
    parser = argparse.ArgumentParser(description="Generate variants from VCF and other genomic data.")
    parser.add_argument("--fai", help="Path to the FAI file.")
    parser.add_argument("--vcf", help="Path to the input VCF file.")
    parser.add_argument("--output", help="Path to the output VCF file.")
    parser.add_argument("--fasta", help="Path to the FASTA file.")
    parser.add_argument("--yaml", help="Path to the YAML file with variant probabilities.")

    args = parser.parse_args()
    main(args.fai, args.vcf, args.output, args.fasta, args.yaml)
