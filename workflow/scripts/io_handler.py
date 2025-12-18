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

    @staticmethod
    def save_subgraph_temp(graph, temp_path: str, sample: str, chromosome: str, max_workers: int = 4):
        """
        Save subgraph to temp file as GFA fragment (S, L, P lines).
        Uses multithreading for larger graphs.
        
        Parameters:
            graph: Graph object to save
            temp_path (str): Path to temp file
            sample (str): Sample name
            chromosome (str): Chromosome identifier
            max_workers (int): Number of threads for parallel processing
        """
        from concurrent.futures import ThreadPoolExecutor, as_completed
        from io import StringIO
        
        node_count = len(graph.nodes)

        # For small graphs, use simple single-threaded approach
        if node_count < 1000 or max_workers <= 1:
            MSpangepopDataHandler._save_subgraph_temp_simple(graph, temp_path, sample, chromosome)
            return

        # Multithreaded approach for larger graphs
        def process_nodes_parallel():
            """Process nodes with chunking"""
            nodes_list = sorted(graph.nodes, key=lambda n: n.id)
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
                chunks.append((len(chunks), chunk_nodes))

            results = [''] * len(chunks)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(process_chunk, chunk): chunk[0] for chunk in chunks}
                for future in as_completed(futures):
                    chunk_idx, chunk_result = future.result()
                    results[chunk_idx] = chunk_result

            return ''.join(results)

        def process_paths_and_edges():
            """Process paths and collect edges"""
            paths_buffer = StringIO()
            unique_edges = set()

            for lineage, path in graph.paths.items():
                if lineage == "Ancestral":
                    continue
                if not path.path_edges:
                    continue

                # Build path string
                path_str = f"{path.path_edges[0].node1.id}+"
                for edge in path.path_edges:
                    orient = "+" if not edge.node2_side else "-"
                    path_str += f",{edge.node2.id}{orient}"

                paths_buffer.write(f"P\t{sample}#{lineage}#chr_{chromosome}\t{path_str}\t*\n")

                # Collect edges
                for edge in path.path_edges:
                    edge_key = (
                        edge.node1.id,
                        "+" if edge.node1_side else "-",
                        edge.node2.id,
                        "+" if not edge.node2_side else "-"
                    )
                    unique_edges.add(edge_key)

            # Convert edges to strings
            edges_buffer = StringIO()
            for n1_id, o1, n2_id, o2 in sorted(unique_edges):
                edges_buffer.write(f"L\t{n1_id}\t{o1}\t{n2_id}\t{o2}\t0M\n")

            return paths_buffer.getvalue(), edges_buffer.getvalue()

        # Execute both tasks in parallel
        with ThreadPoolExecutor(max_workers=2) as executor:
            nodes_future = executor.submit(process_nodes_parallel)
            paths_edges_future = executor.submit(process_paths_and_edges)

            nodes_str = nodes_future.result()
            paths_str, edges_str = paths_edges_future.result()

        # Write to file
        with open(temp_path, 'w') as f:
            f.write(nodes_str)
            f.write(edges_str)
            f.write(paths_str)

    @staticmethod
    def save_subgraph_temp_simple(graph, temp_path: str, sample: str, chromosome: str):
        """
        Simple single-threaded save for small graphs.
        
        Parameters:
            graph: Graph object to save
            temp_path (str): Path to temp file
            sample (str): Sample name
            chromosome (str): Chromosome identifier
        """
        with open(temp_path, 'w') as f:
            # Nodes (S lines)
            for node in sorted(graph.nodes, key=lambda n: n.id):
                f.write(f"S\t{node.id}\t{node}\n")

            # Edges (L lines) - collect unique edges
            unique_edges = set()
            for lineage, path in graph.paths.items():
                if lineage == "Ancestral":
                    continue
                for edge in path.path_edges:
                    edge_key = (
                        edge.node1.id,
                        "+" if edge.node1_side else "-",
                        edge.node2.id,
                        "+" if not edge.node2_side else "-"
                    )
                    unique_edges.add(edge_key)

            for n1_id, o1, n2_id, o2 in sorted(unique_edges):
                f.write(f"L\t{n1_id}\t{o1}\t{n2_id}\t{o2}\t0M\n")

            # Paths (P lines)
            for lineage, path in graph.paths.items():
                if lineage == "Ancestral":
                    continue
                if not path.path_edges:
                    continue

                path_str = f"{path.path_edges[0].node1.id}+"
                for edge in path.path_edges:
                    orient = "+" if not edge.node2_side else "-"
                    path_str += f",{edge.node2.id}{orient}"

                f.write(f"P\t{sample}#{lineage}#chr_{chromosome}\t{path_str}\t*\n")

    @staticmethod
    def merge_temp_files_to_gfa( temp_files: list, output_path: str,sample: str,chromosome: str) -> dict:
        """
        Merge all temp GFA fragments into final GFA.
        Since node IDs are already global, just concatenate and add connecting edges.
        Returns lineage_lengths dict.
        """
        MScompute(f"Merging {len(temp_files)} temp files into final GFA")

        # Collect connecting edges between subgraphs
        connecting_edges = []
        prev_end_id = None

        for tf in temp_files:
            if prev_end_id is not None:
                # Connect previous subgraph's end to this subgraph's start
                connecting_edges.append(f"L\t{prev_end_id}\t+\t{tf.start_node_id}\t+\t0M\n")
            prev_end_id = tf.end_node_id

        # Track path segments for lineage length calculation
        lineage_segments = {}

        with open(output_path, 'w') as out_f:
            # First pass: write all S and L lines, collect P lines
            all_path_lines = []

            for tf in temp_files:
                with open(tf.path, 'r') as in_f:
                    for line in in_f:
                        if line.startswith('S\t') or line.startswith('L\t'):
                            out_f.write(line)
                        elif line.startswith('P\t'):
                            all_path_lines.append(line.strip())

            # Write connecting edges
            for edge_line in connecting_edges:
                out_f.write(edge_line)

            # Merge paths: same lineage across files should be concatenated
            # Parse: P\tsample#lineage#chr_X\tsegments\t*
            paths_by_lineage = {}

            for pline in all_path_lines:
                parts = pline.split('\t')
                path_name = parts[1]  # sample#lineage#chr_X
                segments = parts[2]  # node1+,node2-,...

                # Extract lineage from path name
                lineage = path_name.split('#')[1]

                if lineage not in paths_by_lineage:
                    paths_by_lineage[lineage] = []
                paths_by_lineage[lineage].append(segments)

            # Write merged paths
            for lineage, segment_lists in paths_by_lineage.items():
                # Join segments from consecutive subgraphs
                merged_segments = ",".join(segment_lists)
                path_name = f"{sample}#{lineage}#chr_{chromosome}"
                out_f.write(f"P\t{path_name}\t{merged_segments}\t*\n")

                # Count nodes for lineage length
                lineage_segments[lineage] = merged_segments.count(',') + 1

        MSsuccess(f"Final GFA written: {output_path}")
        return lineage_segments

    @staticmethod
    def write_fasta_from_gfa(gfa_path: str,sample: str,chromosome: str,fasta_folder: str,compress: bool = True):
        """Read the final GFA and write FASTA sequences."""
        import gzip
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        os.makedirs(fasta_folder, exist_ok=True)
        ext = ".fasta.gz" if compress else ".fasta"
        fasta_path = os.path.join(fasta_folder, f"{sample}_chr{chromosome}{ext}")

        MScompute("Building FASTA from GFA")

        # Read node sequences
        node_seqs = {}
        paths = {}  # path_name -> segments string

        with open(gfa_path, 'r') as f:
            for line in f:
                if line.startswith('S\t'):
                    parts = line.strip().split('\t')
                    node_id = int(parts[1])
                    seq = parts[2]
                    node_seqs[node_id] = seq
                elif line.startswith('P\t'):
                    parts = line.strip().split('\t')
                    path_name = parts[1]
                    segments = parts[2]
                    paths[path_name] = segments

        # Reverse complement helper
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

        def rev_comp(s):
            return ''.join(comp.get(b, b) for b in reversed(s))

        # Write FASTA
        open_func = gzip.open if compress else open
        with open_func(fasta_path, 'wt') as out_f:
            for path_name, segments in paths.items():
                seq_parts = []
                for seg in segments.split(','):
                    node_id = int(seg[:-1])
                    orient = seg[-1]
                    node_seq = node_seqs.get(node_id, "")
                    if orient == '-':
                        node_seq = rev_comp(node_seq)
                    seq_parts.append(node_seq)

                full_seq = ''.join(seq_parts)
                record = SeqRecord(Seq(full_seq), id=path_name, description="")
                SeqIO.write(record, out_f, "fasta")

        MSsuccess(f"FASTA written: {fasta_path}")
