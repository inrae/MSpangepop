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
        MScompute(f"Starting to write fasta file with {threads} threads")
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

    @staticmethod
    def save_to_gfav1_1_hybrid(ensemble, file_path, ignore_ancestral=False, max_workers=None):
        """
        Hybrid approach with just-in-time ID assignment during save.
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from io import StringIO
        
        if not ensemble.graphs:
            raise MSerror("No graphs to save")
        
        graph = ensemble.graphs[0]
        
        if max_workers is None:
            max_workers = 4
        
        # ASSIGN IDs 
        MScompute("Assigning node IDs for GFA export...")
        node_to_id = {}
        for idx, node in enumerate(graph.nodes, start=1):
            node.id = idx
            node_to_id[node] = idx
        
        total_nodes = len(graph.nodes)
        total_paths = sum(1 for lineage in graph.paths.keys() 
                        if not (ignore_ancestral and lineage == "Ancestral"))
        MScompute(f"Starting GFA export: {total_nodes} nodes, {total_paths} paths to process")
        
        def process_nodes_parallel():
            """Process nodes with chunking"""
            nodes_list = list(graph.nodes)
            chunk_size = max(1, len(nodes_list) // max_workers)
            
            def process_chunk(chunk_data):
                chunk_idx, chunk = chunk_data
                buffer = StringIO()
                for node in chunk:
                    buffer.write(f"S\t{node.id}\t{node}\n")
                return chunk_idx, buffer.getvalue()
            
            chunks = []
            for i in range(0, len(nodes_list), chunk_size):
                chunk_nodes = nodes_list[i:i + chunk_size]
                chunk_idx = len(chunks)
                chunks.append((chunk_idx, chunk_nodes))
            
            results = [''] * len(chunks)
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_chunk = {executor.submit(process_chunk, chunk_data): chunk_data 
                                for chunk_data in chunks}
                
                for future in as_completed(future_to_chunk):
                    chunk_idx, chunk_result = future.result()
                    results[chunk_idx] = chunk_result
            
            return ''.join(results)
        
        def process_paths_and_edges():
            """Process paths and collect edges"""
            paths_buffer = StringIO()
            unique_edges = set()
            
            for lineage, path in graph.paths.items():
                if ignore_ancestral and lineage == "Ancestral":
                    continue
                
                # Build path representation using node IDs
                path_repr = ""
                if path.path_edges:
                    # Start with first node
                    path_repr = f"{path.path_edges[0].node1.id}+"
                    # Add subsequent nodes
                    for edge in path.path_edges:
                        orientation = '+' if not edge.node2_side else '-'
                        path_repr += f",{edge.node2.id}{orientation}"
                
                paths_buffer.write(f"P\tlineage_{path.lineage}\t{path_repr}\t*\n")
                
                # Collect edges
                for edge in path.path_edges:
                    edge_signature = (
                        edge.node1.id, edge.node1_side,
                        edge.node2.id, edge.node2_side
                    )
                    unique_edges.add(edge_signature)
            
            # Convert edges to strings
            edges_buffer = StringIO()
            for edge_data in unique_edges:
                node1_id, node1_side, node2_id, node2_side = edge_data
                orientation1 = "+" if node1_side else "-"
                orientation2 = "+" if not node2_side else "-"
                edges_buffer.write(f"L\t{node1_id}\t{orientation1}\t{node2_id}\t{orientation2}\t0M\n")
            
            return paths_buffer.getvalue(), edges_buffer.getvalue()
        
        # Execute both tasks in parallel
        MScompute(f"Starting to create buffer with {max_workers} threads")
        
        with ThreadPoolExecutor(max_workers=2) as executor:
            nodes_future = executor.submit(process_nodes_parallel)
            paths_edges_future = executor.submit(process_paths_and_edges)
            
            nodes_str = nodes_future.result()
            paths_str, edges_str = paths_edges_future.result()
        
        # Calculate data sizes
        nodes_size = len(nodes_str)
        edges_size = len(edges_str)
        paths_size = len(paths_str)
        total_size = nodes_size + edges_size + paths_size
        
        MScompute(f"Buffer ready: {total_size:,} characters to save. ({nodes_size:,} S, {edges_size:,} L, {paths_size:,} P)")
        
        # Write to file
        MScompute(f"Writing GFA file...")
        with open(file_path, 'w') as f:
            f.write(nodes_str)
            f.write(edges_str)
            f.write(paths_str)
        
        MSsuccess(f"GFA export completed successfully!")