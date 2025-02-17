import argparse
import itertools
from readfile import read_fasta_gz
from graph_utils import merge_nodes, save_to_gfa
from variants_class import *

class Node:
    """Represents a node in the graph, containing a single nucleotide."""
    def __init__(self, base: str, node_id: int):
        self.id: int = node_id
        self.base: bytes = base.encode("utf-8")  # Store as a single byte to save memory
        self.out_edges: list[Node] = []
        self.in_edges: list[Node] = []

    def connect(self, node: "Node") -> None:
        """Creates a directed edge from this node to another node."""
        self.out_edges.append(node)
        node.in_edges.append(self)

    def __repr__(self) -> str:
        return f"Node(ID={self.id}, {self.base.decode()}, out={len(self.out_edges)}, in={len(self.in_edges)})"

class Graph:
    """Represents a directed graph"""

    _id_counter = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator: itertools.count):
        self.id: int = next(self._id_counter)
        self._node_id_generator: itertools.count = node_id_generator
        self.nodes: list[Node] = []
        self._start_node: Node = None
        self._end_node: Node = None

    def add_node(self, base: str) -> Node:
        """Creates and adds a new node to the graph."""
        node = Node(base, next(self._node_id_generator))
        self.nodes.append(node)
        return node
    
    def __getitem__(self, i: int) -> Node:
        return self.nodes[i]

    def __len__(self) -> int:
        return len(self.nodes)

    def __iadd__(self, other: "Graph") -> "Graph":
        """Concatenates another graph by linking the last node of this one to the first node of the other."""
        if not self.nodes or not other.nodes:
            raise ValueError("Cannot concatenate empty graphs.")
        self.end_node.connect(other.start_node)
        self.nodes.extend(other.nodes)
        self._end_node = other.end_node
        return self

    @property
    def start_node(self) -> Node:
        return self._start_node

    @property
    def end_node(self) -> Node:
        return self._end_node

    def build_from_sequence(self, sequence: str) -> None:
        """Constructs a graph from a nucleotide sequence."""
        if not sequence:
            raise ValueError("❌ MSpangepop -> Cannot build a graph from an empty sequence.")
        
        prev_node = self.add_node(sequence[0])
        self._start_node = prev_node

        for base in sequence[1:]:
            new_node = self.add_node(base)
            prev_node.connect(new_node)
            prev_node = new_node
        
        self._end_node = prev_node

    def add_snp(self) -> None:
        """Adds an SNP at node index 3 if possible."""
        if len(self.nodes) > 3:
            snp = SNP(self.nodes[3])
            new_node = snp.compute_alt_seq(self._node_id_generator)
            self.nodes.append(new_node)

    def __repr__(self) -> str:
        return f"Graph({self.id}, {len(self.nodes)} nodes, ->{self.start_node}-{self.end_node}->)"

class GraphEnsemble:
    """Manages multiple graphs and allows linking them together."""
    
    def __init__(self):
        self.graphs: list[Graph] = []
        self._node_id_generator: itertools.count = itertools.count(1)
    
    def __len__(self) -> int:
        return len(self.graphs)
    
    def __repr__(self) -> str:
        return f"GraphEnsemble({len(self.graphs)} graphs)"
    
    def add_graph(self, graph: Graph) -> None:
        """Adds a graph to the ensemble."""
        self.graphs.append(graph)
    
    def create_empty_graph(self) -> Graph:
        """Creates a new empty graph and adds it to the ensemble."""
        graph = Graph(self._node_id_generator)
        self.add_graph(graph)
        return graph

    @staticmethod
    def concatenate_graphs(ensemble: "GraphEnsemble") -> Graph:
        """Concatenates all graphs in the ensemble and returns the merged graph."""
        if not ensemble.graphs:
            raise ValueError("⚠️ No graphs to concatenate.")

        concatenated_graph = ensemble.graphs[0]
        for graph in ensemble.graphs[1:]:
            concatenated_graph += graph

        merge_nodes(concatenated_graph)
        return concatenated_graph

def main(splited_fasta: str, output_file: str, sample: str, chromosome: str) -> None:
    """Reads a FASTA file and constructs graphs from the sequences."""
    sequences = read_fasta_gz(splited_fasta)
    ensemble = GraphEnsemble()
    
    print(f"\n🔹 MSpangepop -> Constructing sub graphs for {sample}, chr {chromosome}")
    for record in sequences:
        header = record.id  
        seq = str(record.seq)  
        print(f"\t🔹 Handling {header}")
        graph = ensemble.create_empty_graph()
        graph.build_from_sequence(seq)
        graph.add_snp()
    
    print(f"\t✅ MSpangepop -> Constructed {len(sequences)} graphs for {sample}, chr {chromosome}")
    
    print(f"\n🔹 MSpangepop -> Starting to concatenate graphs for {sample}, chr {chromosome}")
    concatenated_graph = GraphEnsemble.concatenate_graphs(ensemble)
    print(f"✅ MSpangepop -> Graph concatenation successful for {sample}, chr {chromosome}")

    print(f"\n🔹 MSpangepop -> Saving graph to {output_file}")
    save_to_gfa(concatenated_graph, output_file)
    print(f"✅ MSpangepop -> Graph saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--sample", required=True, help="Current sample")
    parser.add_argument("--chromosome", required=True, help="Current chromosome")
    
    args = parser.parse_args()
    main(args.splited_fasta, args.output_file, args.sample, args.chromosome)
