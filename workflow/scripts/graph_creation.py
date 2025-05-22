import itertools
import argparse
from bitarray import bitarray
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute
from graph_utils import merge_nodes, save_to_gfa
from variants_lib import *

class Node:
    """Represents a node in the graph."""
    def __init__(self, dna_sequence, node_id: int):
        self.id: int = node_id
        if isinstance(dna_sequence, str):
            self.dna_bases: bitarray = bitarray()
            self.dna_bases.frombytes(dna_sequence.encode("utf-8"))  # Encode string as bitarray
        else:
            self.dna_bases: bitarray = dna_sequence  # Store bitarray directly

        self.outgoing_edges: set["Edge"] = set()  # Store outgoing edges as Edge objects

    @property
    def __decode(self) -> str:
        """Decodes the bitarray back into a DNA sequence string."""
        return self.dna_bases.tobytes().decode("utf-8")

    def __repr__(self) -> str:
        """Returns only the raw DNA sequence as a string."""
        return self.__decode # Using str(Node) will automatically decode the node

    @property
    def reversed(self) -> str:
        string = self.__repr__()
        return string[::-1]

class Edge:
    """Represents a directed edge between two nodes, keeping track of the side of each node involved."""

    def __init__(self, node1: Node, node1_side, node2: Node, node2_side):
        """
        Creates an edge between two nodes.

        Parameters:
        - node1 (Node): First node
        - node1_side (bool): Side of the first node (True for 3', False for 5')
        - node2 (Node): Second node
        - node2_side (bool): Side of the second node (True for 3', False for 5')
        """
        if isinstance(node1_side, str):
            if node1_side == "+" :
                self.node1_side = True
            elif node1_side == "-" :
                self.node1_side = False
        else :
            self.node1_side = node1_side
        if isinstance(node2_side, str):
            if node2_side == "+" :
                self.node2_side = True
            elif node2_side == "-" :
                self.node2_side = False
        else :
            self.node2_side = node2_side

        self.node1 = node1
        self.node2 = node2
        # Add this edge to the outgoing edges of node1
        node1.outgoing_edges.add(self)
    
    def remove(self) -> None:
        """
        Removes this edge from the graph by removing it from node1's outgoing_edges.
        This effectively "detaches" the edge from the graph structure.
        """
        if self in self.node1.outgoing_edges:
            self.node1.outgoing_edges.remove(self)
        else:
            # This should not happen in a well-formed graph, but it's good to handle it
            raise MSerror(f"Edge between node {self.node1.id} and node {self.node2.id} not found in outgoing edges")

class Path:
    """Represents a path in the graph, consisting of a sequence of edges."""

    def __init__(self, lineage, path_edges: list[Edge] = None):
        """
        Initializes a Path.

        Parameters:
        - lineage (int): A unique identifier for the path.
        - path_edges (list[Edge]): A list of Edge objects forming the path.
        """
        self.lineage = lineage
        self.path_edges: list[Edge] = path_edges if path_edges else []
        self.node_count = len(self.path_edges)

    def add_edge(self, edge: Edge) -> None:
        """Adds an edge to the path."""
        self.path_edges.append(edge)
        self.node_count += 1

    def __repr__(self) -> str:
        """Returns a readable representation of the path without repeating nodes."""
        if not self.path_edges:
            return ""  # If there are no edges, return an empty string

        # Start with the first node
        path_repr = f"{self.path_edges[0].node1.id}"

        # Then for each edge, append the node2 id (only if it's not the first node)
        for edge in self.path_edges:
            path_repr += f"{'>' if not edge.node2_side else '<'}{edge.node2.id}"

        return path_repr

    def __str__(self) -> str:
        """Returns a string of the path"""
        if not self.path_edges:
            return ""

        path_repr = f"{self.path_edges[0].node1}"

        for edge in self.path_edges:
            path_repr += f"{edge.node2 if not edge.node2_side else edge.node2.reversed}"

        return path_repr

    def __getitem__(self, i: int) -> Node: # Here we have to reimplement the list functions
        """Returns the node at position i in the path, supporting negative indices."""
        if not self.path_edges:
            raise MSerror("Path is empty")

        if i < 0:
            i = self.node_count + 1 + i  # Convert negative index to positive (like Python lists)

        if i < self.node_count:
            return self.path_edges[i].node1
        elif i == self.node_count:
            return self.path_edges[-1].node2
        else:
            raise MSerror(f"Node index {i} out of bounds for path with {self.node_count + 1} nodes")

    def __iadd__(self, other):
        """Allow concatenation of two paths (or an edge to a path)"""
        if isinstance(other, Edge):
            self.add_edge(other) # If we try to add an edge
        elif isinstance(other, Path):
            for path_edge in other.path_edges:
                self.add_edge(path_edge) # If we try to add a path
        else:
            MSerror("cant += this type to a path")
        return self
    
    def bypass(self, start: int, end: int) -> None:
        """
        Creates a shortcut in the path by creating a direct edge from node at start index
        to node at end index, and removing all edges between them from this path only.
        Does not remove edges from the graph structure to avoid inconsistencies.
        
        Parameters:
        - start (int): Index of the starting node for the bypass
        - end (int): Index of the ending node for the bypass
        """
        if start < 0 or end < 0:
            raise MSerror("Start and end indices must be positive")

        if start > self.node_count or end > self.node_count:
            raise MSerror("Start and end indices out of bounds")
            
        if start == end or start + 1 == end:
            # Nothing to bypass if start and end are the same or adjacent
            return
            
        # Create a new edge from the start node to the end node
        start_node = self[start]
        end_node = self[end]
        
        # Create the new edge - a connection from start to end RETAINING the side of the start node and the end node
        new_edge = Edge(start_node, # Start node
                        self.path_edges[start].node1_side, # Side of the start node
                        end_node, # End node
                        self.path_edges[end-1].node2_side) # Side of the end node
        
        # Replace the set of edges with the new single edge
        # Do NOT remove old edges from graph structure to avoid inconsistencies
        self.path_edges[start:end] = [new_edge]
            
        # Update the node count
        self.node_count = len(self.path_edges)

    def loop(self, start: int, end: int) -> None:
        """
        Creates a loop in the path where the nodes between start and end (inclusive) 
        are traversed twice in sequence.
        
        Examples:
        - loop(1, 2) on A-B-C-D creates A-B-C-B-C-D (nodes B and C are traversed twice)
        - loop(1, 3) on A-B-C-D creates A-B-C-D-B-C-D (nodes B, C, and D are traversed twice)
        - loop(2, 2) on A-B-C-D creates A-B-C-C-D (node C is traversed twice)
        
        Parameters:
        - start (int): Index of the starting node for the loop
        - end (int): Index of the ending node for the loop (inclusive)
        """
        if start < 0 or end < 0:
            raise MSerror("Start and end indices must be positive")
            
        if start > self.node_count or end > self.node_count:
            raise MSerror("Start and end indices out of bounds")
            
        if start > end:
            raise MSerror("Start index must be less than or equal to end index")
        
        # We need to duplicate the edges that connect the nodes in the range [start, end]
        # For nodes at indices [start, end], we need edges at indices [start-1, end-1]
        # But we need to be careful about the boundaries
        
        if start == 0:
            # If we're starting from the first node, we can't duplicate an incoming edge
            # We need to duplicate the edges from start to end
            edges_to_duplicate = self.path_edges[start:end]
        else:
            # Normal case: duplicate edges from (start-1) to (end-1) inclusive
            # This will duplicate the traversal of nodes [start, end]
            edges_to_duplicate = self.path_edges[start-1:end]
        
        # Create copies of the edges to duplicate
        duplicated_edges = []
        for edge in edges_to_duplicate:
            # Create a new edge with the same nodes and sides
            new_edge = Edge(
                edge.node1,      # Same start node
                edge.node1_side, # Same start node side
                edge.node2,      # Same end node  
                edge.node2_side  # Same end node side
            )
            duplicated_edges.append(new_edge)
        
        # Insert the duplicated edges after the original section
        # Insert at position 'end' to place the loop after the original traversal
        self.path_edges[end:end] = duplicated_edges
        
        # Update the node count
        self.node_count = len(self.path_edges)

    def invert(self, start: int, end: int) -> None:
        """
        Create inverted edges in the path between start and end position.
        This function also inverts the side of which the nodes are read between the two positions.
        
        Examples:
        - invert(1, 1) on AB-CD-EF-GH-IJ -> AB-DC-EF-GH-IJ (node CD becomes DC)
        - invert(1, 2) on AB-CD-EF-GH-IJ -> AB-FE-DC-GH-IJ (nodes CD-EF become FE-DC)  
        - invert(1, 3) on AB-CD-EF-GH-IJ -> AB-HG-FE-DC-IJ (nodes CD-EF-GH become HG-FE-DC)
        
        Parameters:
        - start (int): Index of the starting node for the inversion
        - end (int): Index of the ending node for the inversion (inclusive)
        """
        if start < 0 or end < 0: raise MSerror("Start and end indices must be positive")
        if start > self.node_count or end > self.node_count: raise MSerror("Start and end indices out of bounds")
        if start > end: raise MSerror("Start index must be less than or equal to end index")
        
        # Get the nodes in the range [start, end] that need to be inverted
        nodes_to_invert = []
        for i in range(start, end + 1):
            nodes_to_invert.append(self[i])
        
        # Reverse the order of nodes
        nodes_to_invert.reverse()
        
        # Create new edges for the inverted section
        inverted_edges = []
        
        # Handle connection from previous section to inverted section
        if start > 0:
            # Connect from the node before start to the first node of inverted section
            prev_edge = self.path_edges[start - 1]
            first_inverted_node = nodes_to_invert[0]
            
            # Create edge with inverted sides for the first node
            connection_edge = Edge(
                prev_edge.node1,
                prev_edge.node1_side,
                first_inverted_node,
                not prev_edge.node2_side  # Invert the side
            )
            inverted_edges.append(connection_edge)
        
        # Create edges within the inverted section
        for i in range(len(nodes_to_invert) - 1):
            current_node = nodes_to_invert[i]
            next_node = nodes_to_invert[i + 1]
            
            # Find the original edge between these nodes (in reverse order)
            original_edge = None
            for edge in self.path_edges[start:end]:
                if edge.node1 == next_node and edge.node2 == current_node:
                    original_edge = edge
                    break
            
            if original_edge:
                # Create inverted edge (reverse direction and invert sides)
                inverted_edge = Edge(
                    current_node,
                    not original_edge.node2_side,  # Invert the side
                    next_node,
                    not original_edge.node1_side   # Invert the side
                )
            else:
                # Fallback if original edge not found - create with default orientation
                inverted_edge = Edge(
                    current_node,
                    False,  # Default side
                    next_node,
                    True    # Default opposite side
                )
            
            inverted_edges.append(inverted_edge)
        
        # Handle connection from inverted section to next section
        if end < self.node_count:
            # Connect from last inverted node to the next section
            last_inverted_node = nodes_to_invert[-1]
            next_edge = self.path_edges[end]
            
            connection_edge = Edge(
                last_inverted_node,
                not next_edge.node1_side,  # Invert the side
                next_edge.node2,
                next_edge.node2_side
            )
            # This edge will be part of the next section, so we don't add it to inverted_edges
            # Instead, we need to modify the existing edge
        
        # Calculate the range of edges to replace
        edge_start = start - 1 if start > 0 else 0
        edge_end = end if end < self.node_count else self.node_count
        
        # Replace the edges in the inverted section
        if start == 0:
            # If inverting from the beginning, don't include connection from previous
            self.path_edges[0:end] = inverted_edges[1:] if len(inverted_edges) > 1 else []
        else:
            # Replace the section including the connection edge
            self.path_edges[edge_start:edge_end] = inverted_edges
        
        # Update the node count
        self.node_count = len(self.path_edges)
   
class Graph:
    """Represents a directed graph"""

    graph_id_generator = itertools.count(1)  # Unique ID generator for graphs

    def __init__(self, node_id_generator: itertools.count):
        self.id: int = next(self.graph_id_generator)
        self.node_id_generator: itertools.count = node_id_generator
        self.nodes: set[Node] = set()
        self.paths: dict[int, Path] = {}
        self.start_node: Node = None
        self.end_node: Node = None

    def add_new_node(self, dna_sequence: str) -> Node:
        """
        Creates and adds a new node to the graph, return it"""
        node = Node(dna_sequence, next(self.node_id_generator))
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

        # Create connecting edge between self.end_node and other.start_node
        connecting_edge = Edge(self.end_node, True, other.start_node, False)

        # Update end node pointer
        self.end_node = other.end_node

        # Merge paths
        for lineage, other_path in other.paths.items():
            if lineage in self.paths:
                # Get the current path
                current_path = self.paths[lineage]
                if current_path[-1] != connecting_edge.node1:
                    raise MSerror(f"Path discontinuity for lineage {lineage}")

                # Extend the path
                current_path += connecting_edge
                current_path += other_path
            else:
                # If lineage doesn't exist, copy the path and prepend the connecting edge
                copied_path = Path(lineage)
                copied_path += other_path
                self.paths[lineage] = copied_path

        return self

    def build_from_sequence(self, nucleotide_sequence: str, lineages: set) -> None:
        """Constructs a graph and creates the associated paths from a nucleotide sequence."""

        if not nucleotide_sequence:
            raise MSerror("Cannot build a graph from an empty sequence.")

        ancestral_path = Path("Ancestral") # Create the ancestral path

        self.start_node = Node(nucleotide_sequence[0], next(self.node_id_generator)) # Add the first node to the graph
        previous_node = self.start_node
        self.add_node(self.start_node) # Keep the pointer

        for base in nucleotide_sequence[1:]: # For each base
            current_node = self.add_new_node(base) # Create a new node
            ancestral_path += Edge(previous_node, True, current_node, False) # Add an Edge in the ancestral path (+-)
            previous_node = current_node
        self.end_node = previous_node # Keep the pointer on the last node

        self.paths[ancestral_path.lineage] = ancestral_path
        for i in lineages:
            # Create a new path object with a shallow copy of the edges
            copied_path = Path(i, ancestral_path.path_edges.copy())
            self.paths[i] = copied_path

    def __repr__(self) -> str:
        return f"Graph(id={self.id}, Nodes={len(self.nodes)}, Paths={len(self.paths)})"

    @property
    def details(self):
        print(self)
        print("Paths :")
        for path in self.paths.items() :
            print(path)

    def add_del(self, a: int, b: int, affected_lineages: set) -> None:
        """
        Adds a deletion in the specified paths, between positions a and b.
        Creates a bypass in each path to represent the deletion.
        
        Parameters:
        - a (int): Starting position of the deletion
        - b (int): Ending position of the deletion
        - affected_lineages (set): Set of path lineages where the deletion should be applied
    
        """
        if a < 0 or b < 0:
            raise MSerror("Deletion positions must be positive")
            
        if a >= b:
            raise MSerror("End position must be greater than start position")
        
        # Keep track of missing paths
        missing_paths = []
        
        # Apply the deletion (bypass) to each specified path
        for lineage in affected_lineages:
            if lineage in self.paths:
                # Get the path
                path = self.paths[lineage]
                
                # Make sure the path is long enough for the deletion
                if b > path.node_count:
                    raise MSerror(f"Path with lineage {lineage} is too short for deletion between positions {a} and {b}")
                
                # Create a bypass in the path (this will create the deletion)
                path.bypass(a, b)
            else:
                # If a path is not found, add it to missing paths
                missing_paths.append(lineage)
        
        # Warn about missing paths
        if missing_paths:
            print(f"Warning: The following paths were not found in the graph: {missing_paths}")
        
    def add_tdup(self, a: int, b: int, affected_lineages: set) -> None:
        """

        """
        # Make sure the arguments are valid
        if a < 0 or b < 0: raise MSerror("Deletion positions must be positive")
        if a > b: raise MSerror("End must be > start")

        # Keep track of missing paths
        missing_paths = []
        
        # Apply the tdup to each specified path
        for lineage in affected_lineages:
            if lineage in self.paths:
                path = self.paths[lineage]
                if b > path.node_count:
                    raise MSerror(f"Path with lineage {lineage} is too short for deletion between positions {a} and {b}")
                
                # Create a tdup in the path (this will create the deletion)
                path.loop(a, b)
            else:
                missing_paths.append(lineage)
        if missing_paths:
            print(f"Warning: The following paths were not found in the graph: {missing_paths}")

if __name__ == "__main__":
    print("start")
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence("ABCDEFG",[1,2,3])
    graphA.details
    graphA.add_tdup(4, 4, {1, 2})
    graphA.details

    nodeA = Node("|Je|", 1)
    nodeB = Node("|vais|", 2)
    nodeC = Node("|au|", 3)
    nodeD = Node("|concert|", 4)
    nodeE = Node("|ce|", 5)
    nodeF = Node("|soir|", 6)

    path = Path(1, [Edge(nodeA, True, nodeB, False),
                    Edge(nodeB, True, nodeC, False),
                    Edge(nodeC, True, nodeD, False),
                    Edge(nodeD, True, nodeE, False),
                    Edge(nodeE, True, nodeF, False)])
    print(path)
    print(repr(path))
    path.loop(start = 1, end = 3)
    print(path)
    print(repr(path))
    """
    nodeA = Node("|Je|", 1)
    nodeB = Node("|vais|", 2)
    nodeC = Node("|au|", 3)
    nodeD = Node("|concert|", 4)
    nodeE = Node("|ce|", 5)
    nodeF = Node("|soir|", 6)

    path = Path(1, [Edge(nodeA, True, nodeB, False),
                    Edge(nodeB, True, nodeC, True),
                    Edge(nodeC, True, nodeD, False),
                    Edge(nodeD, True, nodeE, True),
                    Edge(nodeE, True, nodeF, False)])
    print(path)
    print(repr(path))
    path.loop(start = 1, end = 4)
    print(path)
    print(repr(path))
    """

'''
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence("ABC",[1,2,3])
    graphA.details
    print()
    graphB = Graph(count)
    graphB.build_from_sequence("DEF",[4])
    graphB.details
    print()
    graphA+=graphB
    graphA.details

class Graph:
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

