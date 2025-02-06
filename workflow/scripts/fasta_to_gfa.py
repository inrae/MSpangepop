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
    
    def __getitem__(self, i):
        return self.nodes[i]

    def __len__(self):
        return len(self.nodes)

    def __add__(self, other):
        """Concatenates two graphs by linking the last node of the first to the first node of the second."""
        if not self.nodes or not other.nodes:
            raise ValueError("Cannot concatenate empty graphs.")
        self.connect_nodes(self.end_node, other.start_node)
        self.nodes.extend(other.nodes)
        self.end_node = other.end_node
        return self

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
    
    def save_to_gfa(self, filename):
        """Saves the concatenated graph to a file in the requested format."""
        if not self.graphs:
            print("⚠️ MSpangepop -> No graphs to save.")
            return

        concatenated_graph = self.graphs[0]
        for graph in self.graphs[1:]:
            concatenated_graph += graph
        
        with open(filename, 'w') as f:
            # Write all nodes
            for node in concatenated_graph.nodes:
                f.write(f"S\t{node.id}\t{node.base.decode()}\n")
            
            # Write all links
            for node in concatenated_graph.nodes:
                for out_node in node.out_edges:
                    f.write(f"L\t{node.id}\t+\t{out_node.id}\t+\t0M\n")
        print(f"✅ MSpangepop -> Graph ensemble saved to {filename}")

    def __repr__(self):
        return f"GraphEnsemble({len(self.graphs)} graphs)"


def main(splited_fasta, output_file):
    """
    Reads a FASTA file with interval sequences, constructs graphs for each, and prints them.
    
    Parameters:
        splited_fasta (str): Path to the FASTA file.
        output_file (str): Path to the output file where graph data will be saved.
    """
    sequences = read_fasta_gz(splited_fasta)
    ensemble = GraphEnsemble()

    for record in sequences:
        header = record.id  
        seq = str(record.seq)  
        print(f"\n🔹 MSpangepop -> Constructing graph for {header}")
        graph = ensemble.create_graph()
        graph.build_from_fasta(seq)
        print("✅ MSpangepop -> ", graph)
        
    # Save the concatenated graph to file
    ensemble.save_to_gfa(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file.")
    
    args = parser.parse_args()
    main(args.splited_fasta, args.output_file)
