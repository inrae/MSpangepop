import itertools
import argparse
from bitarray import bitarray
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute
from graph_utils import merge_nodes, save_to_gfa
from variants_lib import *

class Node:
    """Represents a node in the graph."""

    def __init__(self, base, node_id: int):
        self.id: int = node_id
        if isinstance(base, str):
            self.bases: bitarray = bitarray()
            self.bases.frombytes(base.encode("utf-8"))  # Encode string as bitarray
        else:
            self.base: bitarray = base  # Store bitarray directly
        
        self.outgoing_links: set["Edge"] = set()  # Store outgoing edges as Edge objects

    def decode_base(self) -> str:
        """Decodes the bitarray back into a DNA sequence string."""
        return self.bases.tobytes().decode("utf-8")

    def __repr__(self) -> str:
        """Returns only the raw DNA sequence as a string."""
        return self.decode_base()
    
    @property
    def reversed(self) -> str:
        string = self.__repr__()
        return string[::-1]

class Edge:
    """Represents a directed edge between two nodes, keeping track of the side of each node involved."""
    
    def __init__(self, node1: Node, side1, node2: Node, side2):
        """
        Creates an edge between two nodes.

        Parameters:
        - node1 (Node): First node
        - side1 (bool): Side of the first node (True for 3', False for 5')
        - node2 (Node): Second node
        - side2 (bool): Side of the second node (True for 3', False for 5')
        """
        if isinstance(side1, str):
            if side1 == "+" :
                self.side1 = True 
            elif side1 == "-" :
                self.side1 = False

            if side2 == "+" :
                self.side2 = True 
            elif side2 == "-" :
                self.side2 = False
        else :
            self.side1 = side1
            self.side2 = side2

        self.node1 = node1
        self.node2 = node2

        # Add this edge to the outgoing links of node1
        node1.outgoing_links.add(self)



class Path:
    """Represents a path in the graph, consisting of a sequence of edges."""
    
    def __init__(self, lineage: int, edges: list[Edge] = None):
        """
        Initializes a Path.

        Parameters:
        - lineage (int): A unique identifier for the path.
        - edges (list[Edge]): A list of Edge objects forming the path.
        """
        self.lineage: int = lineage
        self.edges: list[Edge] = edges if edges else []

    def add_edge(self, edge: Edge) -> None:
        """Adds an edge to the path."""
        self.edges.append(edge)

    def __repr__(self) -> str:
        """Returns a readable representation of the path without repeating nodes."""
        if not self.edges:
            return ""  # If there are no edges, return an empty string

        # Start with the first node
        path_repr = f"{self.edges[0].node1.id}"

        # Then for each edge, append the node2 id (only if it's not the first node)
        for edge in self.edges:
            path_repr += f"{'<' if edge.side1 == edge.side2 else '>'}{edge.node2.id}"

        return path_repr
    
    def __str__(self) -> str: 
        """Returns a string of the path"""
        if not self.edges:
            return ""  

        path_repr = f"{self.edges[0].node1}"

        for edge in self.edges:
            path_repr += f"{edge.node2.reversed if edge.side1 == edge.side2 else edge.node2}"

        return path_repr
    
if __name__ == "__main__":
    # Create some nodes
    node_a = Node("Bonjour", node_id=1)
    node_b = Node("Bonjour", node_id=2)
    node_c = Node("Bonjour", node_id=3)
    node_d = Node("Bonjour", node_id=4)

    # Create edges
    edge1 = Edge(node_a, "+", node_b, "-") 
    edge2 = Edge(node_b, True, node_c, True) 
    edge3 = Edge(node_c, True, node_c, False)  
    edge4 = Edge(node_c, False, node_d, True)  

    # Create a path and add edges
    path = Path(lineage=101)
    path.add_edge(edge1)
    path.add_edge(edge2)
    path.add_edge(edge3)
    path.add_edge(edge4)

    # Print the path
    print(repr(path))  
    print(str(path))  

'''
class Graph:
    """Represents a directed graph"""

    _id_counter = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator: itertools.count):
        self.id: int = next(self._id_counter)
        self._node_id_generator: itertools.count = node_id_generator
        self.nodes: set[Node] = []
        self.paths: dict[int, Path] = {} 
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
                raise MSerror(f"Path with lineage {path.lineage} already exists in the graph.")

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
            raise MSerror("Cannot concatenate empty graphs.")

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
                self.paths[lineage] = path #Useless

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
            raise MSerror("Cannot build a graph from an empty sequence.")

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
            raise MSerror("No graphs to concatenate.")

        concatenated_graph = ensemble.graphs[0]
        for graph in ensemble.graphs[1:]:
            concatenated_graph += graph

        return concatenated_graph

def main(splited_fasta: str, output_file: str, sample: str, chromosome: str) -> None:
    """Reads a FASTA file and constructs graphs from the sequences."""
    sequences = MSpangepopDataHandler.read_fasta(splited_fasta)
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
        
        #graph.add_snp(4)  
        #graph.add_deletion(10, 15)
        #graph.add_insertion(20, length= 50)  
        #graph.add_snp(25)
        #graph.add_snp(130)
        #graph.add_inversion(10, 20)
        #graph.add_duplication(30,35)
        #merge_nodes(graph)
            
    MSsuccess(f"Constructed {len(sequences)} graphs for {sample}, chr {chromosome}")
    
    MScompute(f"Starting to concatenate graphs for {sample}, chr {chromosome}")
    concatenated_graph = GraphEnsemble.concatenate_graphs(ensemble)

    MScompute(f"Merging nodes for {sample}, chr {chromosome}")
    merge_nodes(concatenated_graph)

    MScompute(f"Saving graph for {sample}, chr {chromosome}")
    save_to_gfa(concatenated_graph, output_file, sample, chromosome)
    MSsuccess(f"Graph saved for {sample}, chr {chromosome}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Construct a nucleotide graph from an interval FASTA file.")
    parser.add_argument("--splited_fasta", required=True, help="Path to the gzipped FASTA file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file.")
    parser.add_argument("--sample", required=True, help="Current sample")
    parser.add_argument("--chromosome", required=True, help="Current chromosome")
    
    args = parser.parse_args()
    main(args.splited_fasta, args.output_file, args.sample, args.chromosome)'
''' 

