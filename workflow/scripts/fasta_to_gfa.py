import argparse
import itertools
from readfile import read_fasta_gz

class Node:
    """Represents a node in a directed graph, containing a single nucleotide"""
    
    def __init__(self, base: str, node_id):
        self.id = node_id  # Assign a unique ID to the node from GraphEnsemble
        self.base = base.encode("utf-8")  # Store as a single byte to save memory
        self.out_edges = []
        self.in_edges = []

    def add_out_edge(self, node):
        """Creates a directed edge from this node to another node."""
        self.out_edges.append(node)
        node.in_edges.append(self)

    def __repr__(self):
        return f"Node(ID={self.id}, {self.base.decode()}, out={len(self.out_edges)}, in={len(self.in_edges)})"


class Graph:
    """Represents a directed graph where nodes contain single nucleotides (A, T, C, G)."""
    
    _id_counter = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator):
        self.id = next(self._id_counter)  # Assign a unique ID
        self.node_id_generator = node_id_generator  # Shared node ID generator from GraphEnsemble
        self.nodes = []
        self.start_node = None
        self.end_node = None

    def add_node(self, base):
        """Creates a new node and adds it to the graph."""
        node = Node(base, next(self.node_id_generator))
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
        return f"Graph(ID={self.id}, {len(self.nodes)} nodes, start={self.start_node}, end={self.end_node})"


class GraphEnsemble:
    """Manages multiple graphs and allows linking them together."""
    
    def __init__(self):
        self.graphs = []
        self.node_id_generator = itertools.count(1)  # Unique node ID generator across all graphs

    def add_graph(self, graph):
        """Adds a graph to the ensemble."""
        self.graphs.append(graph)

    def create_graph(self):
        """Creates a new graph with a shared node ID generator."""
        graph = Graph(self.node_id_generator)
        self.add_graph(graph)
        return graph

    def link_graphs(self, graph1, graph2):
        """Links two graphs by connecting the last node of graph1 to the first node of graph2."""
        if graph1.end_node and graph2.start_node:
            graph1.connect_nodes(graph1.end_node, graph2.start_node)
            print(f"🔗 MSpangepop -> Linked Graph {graph1.id} to Graph {graph2.id}")
        else:
            print("⚠️ MSpangepop -> Cannot link graphs with missing start/end nodes.")

    def __repr__(self):
        return f"GraphEnsemble({len(self.graphs)} graphs)"


def main(splited_fasta):
    """
    Reads a FASTA file with interval sequences, constructs graphs for each, and prints them.
    
    Parameters:
        splited_fasta (str): Path to the FASTA file.
    """
    sequences = read_fasta_gz(splited_fasta)
    ensemble = GraphEnsemble()

    prev_graph = None
    for record in sequences:
        header = record.id  
        seq = str(record.seq)  
        print(f"\n🔹 MSpangepop -> Constructing graph for {header}")
        graph = ensemble.create_graph()
        graph.build_from_fasta(seq)
        print("✅ MSpangepop -> ", graph)
        
        if prev_graph:
            ensemble.link_graphs(prev_graph, graph)
        prev_graph = graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")

    args = parser.parse_args()
    main(args.splited_fasta)
