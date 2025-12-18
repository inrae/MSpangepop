"""
Author: Lucien Piat
Creation: 20 Oct 2024
Updated: 13 Mar 2025
Institution: INRAe
Project: PangenOak

This script provide three type of classes :
- MSLogger that provide a way to log each operation to stdout wit MSsuccess, MScompute and MSwarning
- MSerror which will raise and error with a custom_traceback
- MSpangepopDataHandler is a class that will handle repeted i/o operations across all scripts like reading or writing json files.

Example : 
> MSsuccess("Hello world")
✅ MSpangepop -> [script_name] Hello world

> sequences = MSpangepopDataHandler.read_fasta(splited_fasta)
"""

import sys
import os
import traceback
import threading

#IDF = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]

IDF = [
    "🌱", "🌿", "🌳", "🌲", "🌻", "🌷", "🪴", "🍄", "🦠", "🧬",
    "🍀", "🌾", "🌼", "🌹", "🪻", "🍁", "🍂", "🌵", "🪹", "🪵"
]

class MSLogger:
    """Base logging class with PID-based identifier."""
    SCRIPT_COL_WIDTH = 20
    THREAD_COL_WIDTH = 10

    def __init__(self, prefix, message):
        self.script = os.path.basename(sys.argv[0][:-3])
        thread_id = threading.get_ident()
        self.identifier = IDF[thread_id % len(IDF)]

        # Fixed-width columns inside the brackets
        script_part = self.script.center(self.SCRIPT_COL_WIDTH)
        thread_part = f"thread:{self.identifier}".center(self.THREAD_COL_WIDTH)

        bracket_text = f"[{script_part}|{thread_part}]"
        print(f"{prefix} {bracket_text} {message}")

class MSsuccess(MSLogger):
    def __init__(self, message):
        super().__init__("✅", message)

class MScompute(MSLogger):
    def __init__(self, message):
        super().__init__("  ", message)

class MSwarning(MSLogger):
    def __init__(self, message):
        super().__init__("⚠️ ", message)

class MSerror(Exception):
    """Custom exception for MSpangepop errors with clean display."""
    def __init__(self, message="An unknown MSpangepop error occurred."):
        self.clean_message = message
        super().__init__(message)  # Store original message

def custom_traceback():
    """Install a custom exception handler for cleaner MSerror display."""
    def handle_mserror(exc_type, exc_value, exc_traceback):
        if exc_type == MSerror:
            print(f"❌ MSpangepop -> Error: {exc_value.clean_message}")
            tb = traceback.extract_tb(exc_traceback)
            if tb:
                last_frame = tb[-1]
                print(f"\t- {os.path.basename(last_frame.filename)}:{last_frame.lineno} in {last_frame.name}()")
        else:
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
    
    sys.excepthook = handle_mserror

custom_traceback()

def get_indent(readable_json):
    if readable_json == True or readable_json == "True":
        return 4
    else :
        return None
    
def process_seed(seed):
    if isinstance(seed, str):
        if seed.lower() == "none":
            return None
        else:
            return min(int.from_bytes(seed.encode(), 'big'), 2**32 - 1)
    elif isinstance(seed, int):
        return max(1, min(seed, 2**32 - 1))
    else:
        raise ValueError("Seed must be either a string or an integer")

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
        import json
        try:
            with open(json_path, 'r') as file:
                return json.load(file)
        except Exception as e:
            raise MSerror(f"Error reading JSON file: {e}")
    
    @staticmethod
    def save_json(data, output_path, readable_json):
        """
        Saves data to a JSON file at the specified output path.
        
        Parameters:
            data (dict): The data to save.
            output_path (str): Path where the JSON file should be saved.
        
        Raises:
            FileReadError: If there is an issue saving the JSON file.
        """
        import json
        try:
            with open(output_path, 'w') as file:
                json.dump(data, file, indent=get_indent(readable_json)) # Set indent to none to reduce json file size
        except Exception as e:
            raise MSerror(f"Error saving JSON file: {e}")
    
    @staticmethod
    def read_fasta(fasta_file):
        """
        Reads a FASTA (possibly gzipped) file and returns a list of sequences.
        
        Parameters:
        fasta_file (str or Path): Path to the FASTA file (compressed or uncompressed).
                                Can be either absolute or relative path.
        
        Returns:
        list: A list of SeqRecord objects from the FASTA file.
        
        Raises:
        FileReadError: If unable to read the FASTA file.
        """
        import gzip
        from Bio import SeqIO
        from pathlib import Path

        # Convert to Path object and resolve to absolute path
        fasta_path = Path(fasta_file).resolve()
        
        # Check if file exists
        if not fasta_path.exists():
            raise MSerror(f"FASTA file not found: {fasta_path}")
        
        if not fasta_path.is_file():
            raise MSerror(f"Path is not a file: {fasta_path}")
        
        # Try reading as compressed file first
        try:
            with gzip.open(fasta_path, "rt") as handle:
                return list(SeqIO.parse(handle, "fasta"))
        except (gzip.BadGzipFile, UnicodeDecodeError):
            # Not a gzip file or corrupted gzip, try uncompressed
            MSwarning("Unable to read as compressed file, trying uncompressed version...")
            try:
                with open(fasta_path, "r") as handle:
                    data = list(SeqIO.parse(handle, "fasta"))
                    MSwarning("Consider compressing the fasta file with bgzip for better performance.")
                    return data
            except Exception as e:
                raise MSerror(f"Unable to read FASTA file '{fasta_path}': {e}")
        except Exception as e:
            # Other errors while reading compressed file
            raise MSerror(f"Error reading FASTA file '{fasta_path}': {e}")
        
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
        import pandas as pd
        try:
            df = pd.read_table(file_path)
            df['cumulative_pb'] = df['pb'].cumsum()
            if not abs(df['cumulative_pb'].iloc[-1] - 1) < 1e-6:
                print(f"⚠️ MSpangepop -> Warning: The cumulative probability of {file_path} is less than 1.")
            return df
        except Exception as e:
            raise MSerror(f"Error reading variant length file {file_path}: {e}")

    @staticmethod
    def read_fai(fai_file):
        """
        Reads a FASTA index (.fai) file and extracts chromosome names and their lengths.

        Parameters:
            fai_file (str): Path to the .fai file.

        Returns:
            numpy.ndarray: A one dimsensional arrray of lenght

        Raises:
            MSerror: If the .fai file does not exist or cannot be read.
        """
        import pandas as pd
        if not os.path.exists(fai_file):
            raise MSerror(f"FAI file not found: {fai_file}")
        try:
            return pd.read_table(fai_file, header=None, usecols=[1], names=["length"])["length"].values
        except Exception as e:
            raise MSerror(f"Error reading FAI file {fai_file}: {e}")      
