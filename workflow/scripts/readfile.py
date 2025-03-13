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
import gzip
import os

def read_json(json_path):
    """
    Reads a JSON file containing chromosome, nodes, edges, and mutations information.
    
    Parameters:
        json_path (str): Path to the JSON file.
        
    Returns:
        dict: A dictionary representation of the JSON data.
    """
    try:
        with open(json_path, 'r') as file:
            data = json.load(file)
        return data
    except Exception as e:
        print(f"❌ MSpangepop -> Error reading JSON file: {e}")
        raise


def save_json(data, output_path):
    """
    Saves the data to a JSON file at the specified output path.
    
    Parameters:
        data (dict): The data to be saved to the JSON file.
        output_path (str): Path where the JSON file should be saved.
    """
    try:
        with open(output_path, 'w') as file:
            json.dump(data, file, indent=4)
    except Exception as e:
        print(f"❌ MSpangepop -> Error saving JSON file: {e}")
        raise


def read_variant_length_file(file_path):
    """
    Reads a variant length distribution file and parses intervals with their probabilities.
    
    Parameters:
        file_path (str): Path to the variant length distribution file.
        
    Returns:
        pandas.DataFrame: DataFrame containing the variant length intervals and cumulative probabilities.
    """
    try:
        df = pd.read_table(file_path)
        df['cumulative_pb'] = df['pb'].cumsum()  # Cumulative sum of probabilities

        # Check if the last cumulative probability is 1
        if not abs(df['cumulative_pb'].iloc[-1] - 1) < 1e-6:
            print(f"⚠️ Warning: The cumulative probability of {file_path}, is less than 1 (this could lead to errors)")
        return df
    except Exception as e:
        print(f"❌ MSpangepop -> Error reading variant length file {file_path}: {e}")
        raise


def read_yaml(yaml_file):
    """
    Reads variant probabilities from a YAML configuration file.
    
    Parameters:
        yaml_file (str): Path to the YAML file containing variant probabilities.
        
    Returns:
        dict: Dictionary of variant types and their associated probabilities.
        
    Raises:
        ValueError: If the sum of the probabilities does not equal 100.
    """
    try:
        with open(yaml_file, 'r') as file:
            variant_probabilities = yaml.safe_load(file)
        
        # Ensure that probabilities sum to 100
        if sum(variant_probabilities.values()) != 100:
            raise ValueError("Sum of variant probabilities in YAML must equal 100.")
        
        return variant_probabilities
    except Exception as e:
        print(f"❌ MSpangepop -> Error reading YAML file: {e}")
        raise


def read_fasta(fasta_file):
    """
    Reads a FASTA (possibly gzipped) file and returns a list of sequences.
    
    Parameters:
        fasta_file (str): Path to the FASTA file (compressed or uncompressed).
        
    Returns:
        list: A list of SeqRecord objects from the FASTA file.
    """
    try:
        with gzip.open(fasta_file, "rt") as handle:
            return list(SeqIO.parse(handle, "fasta"))
    except Exception as e:
        print(f"⚠️ MSpangepop -> Unable to read compressed file, trying uncompressed version...")

        try:
            with open(fasta_file, "r") as handle:
                data = list(SeqIO.parse(handle, "fasta"))
                print("⚠️ MSpangepop -> We recommend compressing the fasta file with bgzip for better performance.")
                return data
        except Exception as e:
            print(f"❌ MSpangepop -> {e}")

        raise IOError("❌ MSpangepop -> Unable to read FASTA file ")
