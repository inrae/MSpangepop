"""
Author: Lucien Piat  
Creation: 4 Feb 2025  
Updated: 4 Feb 2025  
Institution: INRAe  
Project: PangenOak  
"""

import argparse

from readfile import read_fasta_gz

class Node:
    """Represents a node in a directed graph, containing a single nucleotide"""
    
    def __init__(self, base: str):
        self.base = base.encode("utf-8")  # Store as a single byte to save memory
        self.out_edges = []
        self.in_edges = []

    def add_out_edge(self, node):
        """Creates a directed edge from this node to another node."""
        self.out_edges.append(node)
        node.in_edges.append(self)

    def __repr__(self):
        return f"Node({self.base.decode()}, out={len(self.out_edges)}, in={len(self.in_edges)})"


class Graph:
    """Represents a directed graph where nodes contain single nucleotides (A, T, C, G)."""
    
    def __init__(self):
        self.nodes = []
        self.start_node = None
        self.end_node = None

    def add_node(self, base):
        """Creates a new node and adds it to the graph."""
        node = Node(base)
        self.nodes.append(node)
        return node

    def connect_nodes(self, from_node, to_node):
        """Creates a directed edge from `from_node` to `to_node`."""
        from_node.add_out_edge(to_node)

    def build_from_fasta(self, sequence):
        """Constructs a graph from a given nucleotide sequence."""
        if not sequence:
            raise ValueError("❌ MSpangepop -> Cannot build a graph from an empty sequence.")
        
        prev_node = self.add_node(sequence[0])  # First node
        self.start_node = prev_node  # Mark start of graph

        for base in sequence[1:]:
            new_node = self.add_node(base)
            self.connect_nodes(prev_node, new_node)
            prev_node = new_node
        
        self.end_node = prev_node  # Mark end of graph

    def __repr__(self):
        return f"Graph({len(self.nodes)} nodes, start={self.start_node}, end={self.end_node})"


def main(splited_fasta):
    """
    Reads a FASTA file with interval sequences, constructs graphs for each, and prints them.

    Parameters:
        splited_fasta (str): Path to the FASTA file.
    """
    sequences = read_fasta_gz(splited_fasta)

    for record in sequences:
        header = record.id  
        seq = str(record.seq)  
        print(f"\n🔹 MSpangepop -> Constructing graph for {header}")
        graph = Graph()
        graph.build_from_fasta(seq)
        print("✅ MSpangepop -> ", graph)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")

    args = parser.parse_args()
    main(args.splited_fasta)
