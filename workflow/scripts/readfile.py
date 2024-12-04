"""
Author: Lucien Piat
Creation: 20 Oct 2024
Updated: 4 Dec 2024
Institution: INRAe
Project: PangenOak
"""

import yaml
import json
import pandas as pd
from Bio import SeqIO

def read_fai(fai_file):
    """Read a .fai file and return a dictionary of chromosome lengths."""
    try:
        fai_df = pd.read_table(fai_file, header=None, names=["CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"])
        return dict(zip(fai_df["CHROM"], fai_df["LENGTH"]))
    
    except Exception as e:
        print(f"Error reading FAI file: {e}")
        raise

def read_json(json_path):
    """
    Reads a JSON file containing chromosome, nodes, edges, and mutations information.
    
    Parameters:
        json_path (str): Path to the JSON file.
        
    Returns:
        dict: A dictionary representation of the JSON file.
    """
    with open(json_path, 'r') as file:
        data = json.load(file)
    return data

def save_json(data, output_path):
    """
    Saves the data as a JSON file to the specified output path.
    
    Parameters:
        data (dict): The data to save to the file.
        output_path (str): The path where the file should be saved.
    """
    with open(output_path, 'w') as file:
        json.dump(data, file, indent=4)
    print(f"JSON data has been saved to {output_path}")

def read_variant_length_file(file_path):
    """Read length distribution file and parse intervals with probabilities."""
    try:
        df = pd.read_table(file_path)
        df['cumulative_pb'] = df['pb'].cumsum()
        return df
    
    except Exception as e:
        print(f"Error reading variant length file {file_path}: {e}")
        raise

def read_yaml(yaml_file):
    """Reads variant probabilities from a YAML configuration file."""
    try:
        with open(yaml_file, 'r') as file:
            variant_probabilities = yaml.safe_load(file)
        
        # Ensure probabilities sum to 100
        if sum(variant_probabilities.values()) != 100:
            raise ValueError("Sum of variant probabilities in YAML must equal 100.")
        
        return variant_probabilities
    except Exception as e:
        print(f"Error reading YAML file: {e}")
        raise

def read_fasta(input_fasta):
    """Read a FASTA file and return a dictionary of sequences."""
    try:
        return SeqIO.index(input_fasta, "fasta")
    
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        raise

