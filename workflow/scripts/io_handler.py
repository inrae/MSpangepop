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

class MSLogger:
    """Base logging class that prints messages with a standardized prefix."""
    def __init__(self, prefix, message):
        self.script = os.path.basename(sys.argv[0])  # Get the script that is running
        print(f"{prefix} [{self.script}] {message}")

class MSsuccess(MSLogger):
    """Logs a success message."""
    def __init__(self, message):
        super().__init__("✅", message)


class MScompute(MSLogger):
    """Logs a compute-related message."""
    def __init__(self, message):
        super().__init__("🔹", message)

class MSwarning(MSLogger):
    """Logs a warning message."""
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
            # For MSerror, show clean format without full traceback
            print(f"❌ MSpangepop -> Error: {exc_value.clean_message}")
            
            # Optionally show just the relevant line
            tb = traceback.extract_tb(exc_traceback)
            if tb:
                last_frame = tb[-1]
                print(f"\t- {os.path.basename(last_frame.filename)}:{last_frame.lineno} in {last_frame.name}()")
        else:
            # Default handler for other exceptions
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
    
    sys.excepthook = handle_mserror

custom_traceback() # This is used to simplify the traceback message

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
            fasta_file (str): Path to the FASTA file (compressed or uncompressed).
        
        Returns:
            list: A list of SeqRecord objects from the FASTA file.
        
        Raises:
            FileReadError: If unable to read the FASTA file.
        """
        import gzip
        from Bio import SeqIO
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
                raise MSerror(f"Unable to read FASTA file: {e}")
    
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

    @staticmethod
    def write_fasta_multithreaded(graph, sample, chromosome, fasta_folder, threads=4, compress=True):
        """
        Writes lineage sequences from a graph to a FASTA file using multithreading.

        Parameters:
            graph (object): A graph object with a `paths` dictionary of {lineage: Path}.
            sample (str): Sample name for FASTA header.
            chromosome (str): Chromosome identifier.
            fasta_folder (str): Directory to write the output FASTA file.
            threads (int): Number of worker threads to use for record generation.
            compress (bool): Whether to gzip compress the output file. Default: True
        """
        import gzip
        from Bio import SeqIO
        from concurrent.futures import ThreadPoolExecutor
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        os.makedirs(fasta_folder, exist_ok=True)
        ext = ".fasta.gz" if compress else ".fasta"
        fasta_output_path = os.path.join(fasta_folder, f"{sample}_chr{chromosome}{ext}")

        def generate_seqrecord(lineage, path):
            if lineage == "Ancestral":
                return None  # Skip "Ancestral" lineage
            try:
                sequence_str = str(path)
                record_id = f"sample_{sample}#lineage_{lineage}#chr_{chromosome}"
                return SeqRecord(Seq(sequence_str), id=record_id, description="")
            except Exception as e:
                MSwarning(f"Failed to convert path for lineage '{lineage}': {e}")
                return None

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [
                executor.submit(generate_seqrecord, lineage, path)
                for lineage, path in graph.paths.items()
                if lineage != "Ancestral"
            ]

            # Open gzipped or plain text output file
            open_func = gzip.open if compress else open
            with open_func(fasta_output_path, "wt") as out_handle:
                count = 0
                for future in futures:
                    record = future.result()
                    if record:
                        SeqIO.write(record, out_handle, "fasta")
                        count += 1

        MSsuccess(f"Wrote {count} lineage FASTA sequences to {fasta_output_path}")