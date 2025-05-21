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

    @property
    def __decode(self) -> str:
        """Decodes the bitarray back into a DNA sequence string."""
        return self.bases.tobytes().decode("utf-8")

    def __repr__(self) -> str:
        """Returns only the raw DNA sequence as a string."""
        return self.__decode # Using str(Node) will automaticly decode the node
    
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
        else :
            self.side1 = side1
        if isinstance(side2, str):
            if side2 == "+" :
                self.side2 = True 
            elif side2 == "-" :
                self.side2 = False
        else :
            self.side2 = side2

        self.node1 = node1
        self.node2 = node2
        # Add this edge to the outgoing links of node1
        node1.outgoing_links.add(self)

class Path:
    """Represents a path in the graph, consisting of a sequence of edges."""
    
    def __init__(self, lineage, edges: list[Edge] = None):
        """
        Initializes a Path.

        Parameters:
        - lineage (int): A unique identifier for the path.
        - edges (list[Edge]): A list of Edge objects forming the path.
        """
        self.lineage = lineage
        self.edges: list[Edge] = edges if edges else []
        self.nodes = len(self.edges)

    def add_edge(self, edge: Edge) -> None:
        """Adds an edge to the path."""
        self.edges.append(edge)
        self.nodes += 1

    def __repr__(self) -> str:
        """Returns a readable representation of the path without repeating nodes."""
        if not self.edges:
            return ""  # If there are no edges, return an empty string

        # Start with the first node
        path_repr = f"{self.edges[0].node1.id}"

        # Then for each edge, append the node2 id (only if it's not the first node)
        for edge in self.edges:
            path_repr += f"{'>' if not edge.side2  else '<'}{edge.node2.id}"

        return path_repr
    
    def __str__(self) -> str: 
        """Returns a string of the path"""
        if not self.edges:
            return ""  

        path_repr = f"{self.edges[0].node1}"

        for edge in self.edges:
            path_repr += f"{edge.node2 if not edge.side2 else edge.node2.reversed}"

        return path_repr
    
    def __getitem__(self, i: int) -> Node: # Ici on est obligé de réimplémanter les fonction des listes
        """Returns the node at position i in the path, supporting negative indices."""
        if not self.edges:
            raise MSerror("Path is empty")

        if i < 0:
            i = self.nodes + 1 + i  # Convert negative index to positive (like Python lists)

        if i < self.nodes:
            return self.edges[i].node1
        elif i == self.nodes:
            return self.edges[-1].node2
        else:
            raise MSerror(f"Node index {i} out of bounds for path with {self.nodes + 1} nodes")

    
    def __iadd__(self, other):
        """Allow concatenation of two paths (or a edge to a path)"""
        if isinstance(other, Edge):
            self.add_edge(other) # Si on essayer d'ajouter un edge
        elif isinstance(other, Path):
            for edge_ in other.edges:
                self.add_edge(edge_) # Si on essaye d'ajouter un path
        else:
            MSerror("cant += this type to a path")
        return self


### YOU STOPED HERE

class Graph:
    """Represents a directed graph"""

    _id_counter = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator: itertools.count):
        self.id: int = next(self._id_counter)
        self._node_id_generator: itertools.count = node_id_generator
        self.nodes: set[Node] = set()
        self.paths: dict[int, Path] = {} 
        self._start_node: Node = None
        self._end_node: Node = None

    def add_new_node(self, base: str) -> Node:
        """
        Creates and adds a new node to the graph, return it"""
        node = Node(base, next(self._node_id_generator))
        self.nodes.add(node)
        return node
    
    def add_node(self, node: Node) -> None:
        """Adds a new node to the graph."""
        self.nodes.add(node)

    def __iadd__(self, other: "Graph") -> "Graph":
        """Merges another graph into this one."""

        if not self.nodes or not other.nodes:
            raise MSerror("Cannot concatenate empty graphs.")

        # Merge all nodes
        self.nodes.update(other.nodes)

        # Create connecting edge between self._end_node and other._start_node
        connecting_edge = Edge(self._end_node, True, other._start_node, False)

        # Update end node pointer
        self._end_node = other._end_node

        # Merge paths
        for lineage, o_path in other.paths.items():
            if lineage in self.paths:
                # Get the current path
                current_path = self.paths[lineage]
                if current_path[-1] != connecting_edge.node1:
                    raise MSerror(f"Path discontinuity for lineage {lineage}")

                # Extend the path
                current_path += connecting_edge
                current_path += o_path
            else:
                # If lineage doesn't exist, copy the path and prepend the connecting edge
                new_path = Path(lineage)
                new_path += o_path
                self.paths[lineage] = new_path

        return self

    def build_from_sequence(self, sequence: str, lineages: set) -> None:
        """Constructs a graph and creates the associated paths from a nucleotide sequence."""

        if not sequence:
            raise MSerror("Cannot build a graph from an empty sequence.")

        ancestral_path = Path("Ancestral") # On crée le chemin ancestral

        self._start_node = Node(sequence[0], next(self._node_id_generator)) # On ajoute la premiére node au grph
        prev_node = self._start_node
        self.add_node(self._start_node) # On garde le pointeur

        for base in sequence[1:]: # Pour chaque bases
            new_node = self.add_new_node(base) # On crée une nouvelle node
            ancestral_path += Edge(prev_node, True, new_node, False) # On ajoute un Edge dans le chemin ancestral (+-)
            prev_node = new_node 
        self._end_node = prev_node # On garde le pointeur sur la dernière node

        self.paths[ancestral_path.lineage] = ancestral_path
        for i in lineages:
            # Create a new path object with a shallow copy of the edges
            path_copy = Path(i, ancestral_path.edges.copy())
            self.paths[i] = path_copy

    def __repr__(self) -> str:
        return f"Graph(id={self.id}, Nodes={len(self.nodes)}, Paths={len(self.paths)})"
    
    @property
    def details(self):
        print(self)
        print("Paths :")
        for path in self.paths.items() :
            print(path)


if __name__ == "__main__":
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence(sequence = "ABC", lineages = [1,2,3])
    graphA.details
    print()
    graphB = Graph(count)
    graphB.build_from_sequence(sequence = "DEF", lineages = [4])
    graphB.details
    print()
    graphA+=graphB
    graphA.details


'''
class Graph:

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

