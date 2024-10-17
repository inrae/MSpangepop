from bed2vcf import read_fa, get_seq
from randombed import generate_type
import argparse
import pandas as pd
from multiprocessing import Process

def sv_vcf(vcf_file, fasta_file, fai_file, yaml_file, output_name):
    """
    Generate a VCF file with structural variants based on an input VCF and configuration YAML.
    """
    # Read the input VCF file into a DataFrame
    vcf_df = pd.read_table(vcf_file, sep="\t")

    # Get the number of variants in the VCF file
    num_variants = len(vcf_df)

    # Generate a BED DataFrame of structural variants based on the configuration
    bed_df = generate_type(num_variants, yaml_file, fai_file)

    # Read the FASTA file for sequence extraction
    reference_sequence = read_fa(fasta_file)

    # Extract sequences for the variants and create the final VCF file
    get_seq(vcf_df, bed_df, reference_sequence, output_name)

# Argument parser to handle command-line inputs
parser = argparse.ArgumentParser(description='Create a final VCF file with structural variants.')
parser.add_argument('-v', '--vcf', type=str, required=True, help='Input VCF file generated with msprime.')
parser.add_argument('-fa', '--fasta', type=str, required=True, help='Reference FASTA file for the variants.')
parser.add_argument('-fai', '--fai', type=str, required=True, help='Samtools index of the reference FASTA.')
parser.add_argument('-y', '--yaml', type=str, required=True, help='YAML configuration for structural variants.')
parser.add_argument('-o', '--outName', type=str, required=True, help='Output name for the final VCF file.')

if __name__ == '__main__':
    # Parse command-line arguments
    args = parser.parse_args()
    args_list = list(vars(args).values())

    # Create a separate process to handle the variant generation
    process = Process(target=sv_vcf, args=args_list)
    process.start()
    process.join()
