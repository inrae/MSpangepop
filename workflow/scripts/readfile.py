"""
Author: Lucien Piat
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak
"""

import yaml
import pandas as pd
from Bio import SeqIO


def read_yaml(file_path):
    """Read a YAML file and ensure that percentages sum to 100."""
    try:
        with open(file_path, 'r') as stream:
            data = yaml.safe_load(stream)
        
        if sum(data.values()) != 100:
            raise ValueError("Sum of values in the config file MUST be 100.")
        
        return data
    
    except Exception as e:
        print(f"Error reading YAML file: {e}")
        raise

def read_vcf(vcf_file):
    """Read a VCF file, format it with custom column names, and extract chromosome, position, and sample data."""
    try:
        # Read the first line to determine the number of columns
        with open(vcf_file, 'r') as file:
            first_line = file.readline().strip()
        
        # Split the first line to get the number of columns
        cols = first_line.split('\t')
        num_samples = len(cols) - 9  # Subtract the first 9 standard VCF columns
        
        # Create a dynamic list of sample column names
        sample_cols = [f"SAMPLE{i + 1}" for i in range(num_samples)]
        all_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_cols
        
        # Read the entire VCF file with the appropriate column names
        vcf_df = pd.read_table(vcf_file, sep="\t", header=None, names=all_cols)
        vcf_df['CHROM'] = vcf_df['CHROM'].astype(str)
        vcf_df['POS'] = vcf_df['POS'].astype(int)

        return vcf_df[['CHROM', 'POS'] + sample_cols]
    
    except Exception as e:
        print(f"Error reading VCF file: {e}")
        raise

def read_fai(fai_file):
    """Read a .fai file and return a dictionary of chromosome lengths."""
    try:
        fai_df = pd.read_table(fai_file, header=None, names=["CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"])
        return dict(zip(fai_df["CHROM"], fai_df["LENGTH"]))
    
    except Exception as e:
        print(f"Error reading FAI file: {e}")
        raise

def read_variant_length_file(file_path):
    """Read length distribution file and parse intervals with probabilities."""
    try:
        df = pd.read_table(file_path)
        df['cumulative_pb'] = df['pb'].cumsum()
        return df
    
    except Exception as e:
        print(f"Error reading variant length file {file_path}: {e}")
        raise

def read_fasta(input_fasta):
    """Read a FASTA file and return a dictionary of sequences."""
    try:
        return SeqIO.index(input_fasta, "fasta")
    
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        raise
