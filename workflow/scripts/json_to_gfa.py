import itertools
import argparse
from readfile import read_fasta
from graph_utils import merge_nodes, save_to_gfa
from variants_lib import *

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

class Path:
    """Represents a path in a graph, consisting of a sequence of nodes."""
    def __init__(self, lineage: str, ancesters: set, nodes: list[Node] = []):
        self.lineage: int = lineage
        self.nodes: list[Node] = nodes  # List of nodes in the path
        self.ancesters: set = ancesters  # A set representing the lineage of the path

    def __repr__(self) -> str:
        """String representation of the path."""
        node_ids = [node.id for node in self.nodes]
        return f"Path(path={self.name}, nodes={node_ids})"

class Graph:
    """Represents a directed graph"""

    _id_counter = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator: itertools.count):
        self.id: int = next(self._id_counter)
        self._node_id_generator: itertools.count = node_id_generator
        self.nodes: list[Node] = []
        self.paths: dict[int, Path] = {}  # Dictionary to store paths by lineage
        self._start_node: Node = None
        self._end_node: Node = None

    def add_node(self, base: str) -> Node:
        """Creates and adds a new node to the graph."""
        node = Node(base, next(self._node_id_generator))
        self.nodes.append(node)
        return node

    def get_path(self, lineage: int) -> Path:
        """Retrieves a path by its lineage."""
        return self.paths.get(lineage, None)
    
    def initialize_paths(self, paths: set[Path]) -> None:
        """
        Initializes the given set of paths, ensuring that:
        - Each path receives an ordered list of nodes from the graph.
        - Each path has a unique lineage in the graph.
        """
        for path in paths:
            if path.lineage in self.paths:
                raise ValueError(f"Path with lineage {path.lineage} already exists in the graph.")

            # Assign the ordered nodes from the graph to the path
            path.nodes = self.nodes[:]  # Shallow copy of the current nodes list
            self.paths[path.lineage] = path

    def __getitem__(self, i: int) -> Node:
        return self.nodes[i]

    def __len__(self) -> int:
        return len(self.nodes)

    def __iadd__(self, other: "Graph") -> "Graph":
        """Concatenates another graph by linking the last node of this one to the first node of the other.
        Also merges paths that share the same lineage.
        """
        if not self.nodes or not other.nodes:
            raise ValueError("Cannot concatenate empty graphs.")

        # Connect the end node of this graph to the start node of the other graph
        self.end_node.connect(other.start_node)
        self.nodes.extend(other.nodes)
        self._end_node = other.end_node

        # Merge paths based on lineage
        for lineage, path in other.paths.items():
            if lineage in self.paths:
                # Concatenate paths with the same lineage
                self.paths[lineage].nodes.extend(path.nodes)
                self.paths[lineage].ancesters.update(path.ancesters)
            else:
                # Add new paths that don't exist in self
                self.paths[lineage] = path

        return self

    def __repr__(self) -> str:
        return f"Graph({self.id}, {len(self.nodes)} nodes, {len(self.paths)} paths, ->{self.start_node}-{self.end_node}->)"

    @property
    def start_node(self) -> Node:
        return self._start_node

    @property
    def end_node(self) -> Node:
        return self._end_node

    def build_from_sequence(self, sequence: str) -> None:
        """Constructs a graph and creates an associated path from a nucleotide sequence."""
        if not sequence:
            raise ValueError("❌ MSpangepop -> Cannot build a graph from an empty sequence.")

        prev_node = self.add_node(sequence[0])
        self._start_node = prev_node

        for base in sequence[1:]:
            new_node = self.add_node(base)
            prev_node.connect(new_node)
            prev_node = new_node

        self._end_node = prev_node
    
    def add_snp(self, idx: int) -> None:
        snp = SNP(A=self.nodes[idx-1], B=self.nodes[idx], D=self.nodes[idx+1])
        new_node = snp.compute_alt_seq(self._node_id_generator)
        self.nodes.append(new_node)

    def add_deletion(self, start_idx: int, end_idx: int) -> None:
        deletion = DEL(A=self.nodes[start_idx], D=self.nodes[end_idx])
        deletion.compute_alt_seq()

    def add_insertion(self, idx: int, length) -> None:
        insertion = INS(A=self.nodes[idx], D=self.nodes[idx + 1], length=length)
        inserted_nodes = insertion.compute_alt_seq(self._node_id_generator)
        self.nodes = self.nodes + inserted_nodes 
    
    def add_inversion(self, start_idx: int, end_idx: int) -> None:
        inversion = INV(A=self.nodes[start_idx], B=self.nodes[start_idx+1], C=self.nodes[end_idx-1], D=self.nodes[end_idx])
        inversion.compute_alt_seq()

    def add_duplication(self, start_idx: int, end_idx: int) -> None:
        duplication = DUP(A=self.nodes[start_idx], D=self.nodes[end_idx])
        duplication.compute_alt_seq()

        
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

        return concatenated_graph


def main(splited_fasta: str, output_file: str, sample: str, chromosome: str) -> None:
    """Reads a FASTA file and constructs graphs from the sequences."""
    sequences = read_fasta(splited_fasta)
    ensemble = GraphEnsemble()
    
    print(f"\n🔹 MSpangepop -> Constructing sub graphs for {sample}, chr {chromosome}")
    for record in sequences:
        header = record.id  
        seq = str(record.seq)  
        print(f"\t🔹 Handling {header}")
        graph = ensemble.create_empty_graph()
        graph.build_from_sequence(seq)

        path = Path(lineage=1, ancesters={1, 2, 3})
        graph.initialize_paths({path})
        
        graph.add_snp(3)  # Add SNP at index 3
        graph.add_deletion(5, 10)  # Delete bases between 5 and 10
        graph.add_insertion(15, 20)  # Insert 10 random bases at position 15
        graph.add_insertion(110, 5)
        graph.add_inversion(40, 80)
        graph.add_duplication(25,35)
        graph.add_snp(10)
        graph.add_snp(5)
        graph.add_snp(8)
        graph.add_snp(50)
        graph.add_snp(83)


    print(f"\t✅ MSpangepop -> Constructed {len(sequences)} graphs for {sample}, chr {chromosome}")
    
    print(f"\n🔹 MSpangepop -> Starting to concatenate graphs for {sample}, chr {chromosome}")
    concatenated_graph = GraphEnsemble.concatenate_graphs(ensemble)
    print(f"✅ MSpangepop -> Graph concatenation  for {sample}, chr {chromosome}")

    print(f"\n🔹 MSpangepop -> Merge nodes for {sample}, chr {chromosome}")
    merge_nodes(concatenated_graph)
    print(f"✅ MSpangepop -> Nodes merged successfuly for {sample}, chr {chromosome}")

    print(f"\n🔹 MSpangepop -> Saving graph to {output_file}")
    save_to_gfa(concatenated_graph, output_file, sample, chromosome)
    print(f"✅ MSpangepop -> Graph saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--sample", required=True, help="Current sample")
    parser.add_argument("--chromosome", required=True, help="Current chromosome")
    
    args = parser.parse_args()
    main(args.splited_fasta, args.output_file, args.sample, args.chromosome)