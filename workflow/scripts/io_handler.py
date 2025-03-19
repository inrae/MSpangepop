"""
Author: Lucien Piat
Creation: 20 Oct 2024
Updated: 13 Mar 2025
Institution: INRAe
Project: PangenOak
"""

import json
import yaml
import pandas as pd
from Bio import SeqIO
import gzip
import sys
import os

import sys
import os

class MSLogger:
    """Base logging class that prints messages with a standardized prefix."""
    def __init__(self, prefix, message):
        self.script = os.path.basename(sys.argv[0])  # Get the script that is running
        print(f"{prefix} [{self.script}] {message}")

class MSsuccess(MSLogger):
    """Logs a success message."""
    def __init__(self, message):
        super().__init__("✅ MSpangepop ->", message)


class MScompute(MSLogger):
    """Logs a compute-related message."""
    def __init__(self, message):
        super().__init__("🔹MSpangepop ->", message)


class MSwarning(MSLogger):
    """Logs a warning message."""
    def __init__(self, message):
        super().__init__("⚠️ MSpangepop ->", message)

class MSerror(Exception):
    """Custom exception for MSpangepop errors."""
    def __init__(self, message):
        self.script = os.path.basename(sys.argv[0])  # Get the script that triggered the error
        self.message = f"❌ MSpangepop -> [{self.script}] {message}"
        super().__init__(self.message)

class MSpangepopDataHandler:
    """
    Handles file operations and processing of variant length distributions and probabilities.
    """
    
    @staticmethod
    def read_json(json_path):
        """
        Reads a JSON file and returns its contents as a dictionary.
        
        Parameters:
            json_path (str): Path to the JSON file.
        
        Returns:
            dict: Dictionary representation of the JSON data.
        
        Raises:
            FileReadError: If there is an issue reading the JSON file.
        """
        try:
            with open(json_path, 'r') as file:
                return json.load(file)
        except Exception as e:
            MSerror(f"Error reading JSON file: {e}")
    
    @staticmethod
    def save_json(data, output_path, readble_json = False):
        """
        Saves data to a JSON file at the specified output path.
        
        Parameters:
            data (dict): The data to save.
            output_path (str): Path where the JSON file should be saved.
        
        Raises:
            FileReadError: If there is an issue saving the JSON file.
        """
        if readble_json:
            indent = 4
        else :
            indent = None
        try:
            with open(output_path, 'w') as file:
                json.dump(data, file, indent=indent) # Set indent to none to reduce json file size
        except Exception as e:
            MSerror(f"Error saving JSON file: {e}")
    
    @staticmethod
    def read_yaml(yaml_file):
        """
        Reads variant probabilities from a YAML configuration file.
        
        Parameters:
            yaml_file (str): Path to the YAML file containing variant probabilities.
        
        Returns:
            dict: Dictionary of variant types and their associated probabilities.
        
        Raises:
            ValueError: If the sum of the probabilities does not equal 100.
            FileReadError: If there is an issue reading the YAML file.
        """
        try:
            with open(yaml_file, 'r') as file:
                variant_probabilities = yaml.safe_load(file)
            if sum(variant_probabilities.values()) != 100:
                raise ValueError("Sum of variant probabilities in YAML must equal 100.")
            return variant_probabilities
        except Exception as e:
            MSerror(f"Error reading YAML file: {e}")
    
    @staticmethod
    def read_fasta(fasta_file):
        """
        Reads a FASTA (possibly gzipped) file and returns a list of sequences.
        
        Parameters:
            fasta_file (str): Path to the FASTA file (compressed or uncompressed).
        
        Returns:
            list: A list of SeqRecord objects from the FASTA file.
        
        Raises:
            FileReadError: If unable to read the FASTA file.
        """
        try:
            with gzip.open(fasta_file, "rt") as handle:
                return list(SeqIO.parse(handle, "fasta"))
        except Exception:
            MSwarning("Unable to read compressed file, trying uncompressed version...")
            try:
                with open(fasta_file, "r") as handle:
                    data = list(SeqIO.parse(handle, "fasta"))
                    MSwarning("Consider compressing the fasta file with bgzip for better performance.")
                    return data
            except Exception as e:
                MSerror(f"Unable to read FASTA file: {e}")
    
    @staticmethod
    def read_variant_length_file(file_path):
        """
        Reads a variant length distribution file and parses intervals with their probabilities.
        
        Parameters:
            file_path (str): Path to the variant length distribution file.
        
        Returns:
            pandas.DataFrame: DataFrame containing the variant length intervals and cumulative probabilities.
        
        Raises:
            FileReadError: If there is an issue reading the file.
        """
        try:
            df = pd.read_table(file_path)
            df['cumulative_pb'] = df['pb'].cumsum()
            if not abs(df['cumulative_pb'].iloc[-1] - 1) < 1e-6:
                print(f"⚠️ MSpangepop -> Warning: The cumulative probability of {file_path} is less than 1.")
            return df
        except Exception as e:
            MSerror(f"Error reading variant length file {file_path}: {e}")
