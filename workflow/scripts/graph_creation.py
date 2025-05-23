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
        
    def __repr__(self):
        if self.node1_side: 
            side1 = "+"
        else:
            side1 = "-"
        if self.node2_side: 
            side2 = "+"
        else:
            side2 = "-"
        return f"({self.node1.id}{side1} -> {self.node2.id}{side2})"

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
        Inverts a section of the path by reversing node order and flipping edge orientations.
        This represents a genomic inversion where DNA is read in reverse.
        
        For path with nodes [0,1,2,3,4] and edges [0>1, 1>2, 2>3, 3>4]:
        - invert(1,3) reverses nodes [1,2,3] to [3,2,1] and flips all affected edges
        - Result: [0<3, 3>2, 2<1, 1>4] (nodes 1,2,3 are now read backwards)
        
        Parameters:
        - start (int): Index of first node to invert (1-based for edges, 0-based for nodes)
        - end (int): Index of last node to invert (inclusive)
        """
        # Input validation
        if start <= 0 or end <= 0: 
            raise MSerror("Start and end indices must be positive")
        if start > self.node_count or end + 1 > self.node_count: 
            raise MSerror("Start and end indices out of bounds")
        if start > end: 
            raise MSerror("Start index must be less than or equal to end index")
        
        # SINGLE NODE INVERSION
        if start == end:
            # For single node inversion, we flip the orientations of edges connecting to this node
            
            # Flip the edge leading TO this node (if it exists)
            if start > 0:
                # path_edges[start-1] is the edge that leads TO node at index 'start'
                # We use start-1 because edge indices are offset by 1 from node indices
                edge_to_node = self.path_edges[start - 1]
                new_edge_to_node = Edge(
                    edge_to_node.node1,
                    edge_to_node.node1_side,        # Previous node exit side unchanged
                    edge_to_node.node2,             # Target node unchanged
                    not edge_to_node.node2_side     # FLIP: invert how target node is entered
                )
                self.path_edges[start - 1] = new_edge_to_node
            
            # Flip the edge leading FROM this node (if it exists)
            if start < self.node_count:
                # path_edges[start] is the edge that leads FROM node at index 'start'
                edge_from_node = self.path_edges[start]
                new_edge_from = Edge(
                    edge_from_node.node1,
                    not edge_from_node.node1_side,  # FLIP: invert how source node is exited
                    edge_from_node.node2,
                    edge_from_node.node2_side       # Next node entry side unchanged
                )
                self.path_edges[start] = new_edge_from
            return
        
        # MULTI-NODE INVERSION
        
        # Step 1: Identify boundary edges
        # These are the edges that connect the inversion to the rest of the path
        first_edge = self.path_edges[start - 1]  # Edge leading TO first inverted node
        last_edge = self.path_edges[end]         # Edge leading FROM last inverted node
        
        # Step 2: Extract edges within the inversion
        # path_edges[start:end] gives us edges BETWEEN the inverted nodes
        # We use start:end (not start-1:end+1) because we want internal edges only
        original_inversion = self.path_edges[start:end]
        
        # Step 3: Reverse and flip the internal edges
        original_inversion.reverse()  # Reverse order for inversion
        inverted_inversion = []
        
        for edge in original_inversion:
            # Completely flip the edge: swap nodes and keep their sides
            # This reverses the direction while maintaining proper orientation
            flipped_edge = Edge(
                edge.node2,      # Destination becomes source
                edge.node2_side, # Keep the side orientation
                edge.node1,      # Source becomes destination  
                edge.node1_side  # Keep the side orientation
            )
            inverted_inversion.append(flipped_edge)
        
        # Step 4: Create new boundary edges
        # These connect the inverted section to the unchanged parts of the path
        
        # New edge connecting previous section to first inverted node
        new_first_edge = Edge(
            first_edge.node1,                      # Source: unchanged
            first_edge.node1_side,                 # Source side: unchanged
            inverted_inversion[0].node1,          # Target: first node of inverted section
            not inverted_inversion[0].node1_side   # Target side: FLIPPED to invert reading
        )
        
        # New edge connecting last inverted node to next section
        new_last_edge = Edge(
            inverted_inversion[-1].node2,          # Source: last node of inverted section
            not inverted_inversion[-1].node2_side, # Source side: FLIPPED to invert reading
            last_edge.node2,                       # Target: unchanged
            last_edge.node2_side                   # Target side: unchanged
        )
        
        # Step 5: Replace the entire affected section
        # We replace edges from (start-1) to (end+1) exclusive, which is (start-1) to (end) inclusive
        # This replaces: [boundary_in, internal_edges..., boundary_out]
        final_edges = [new_first_edge] + inverted_inversion + [new_last_edge]
        self.path_edges[start - 1:end + 1] = final_edges
        
        # Update node count (should remain the same for inversions)
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

    def add_del(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds a deletion in the specified paths, between positions a and b.
        Creates a bypass in each path to represent the deletion.
        
        Parameters:
        - a (int): Starting position of the deletion
        - b (int): Ending position of the deletion  
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        # Validate input parameters
        if a < 0 or b < 0: raise MSerror("Deletion positions must be positive")
        if a >= b: raise MSerror("End position must be greater than start position")
        
        # Handle single lineage or collection
        lineages_to_process = [affected_lineages] if isinstance(affected_lineages, (int, str)) else affected_lineages
        missing_paths = []
        
        # Apply deletion to each specified path
        for lineage in lineages_to_process:
            if lineage in self.paths:
                path = self.paths[lineage]
                # Check path length before operation
                if b > path.node_count:
                    raise MSerror(f"Path {lineage} too short for deletion at positions {a}-{b}")
                # Apply bypass operation (creates deletion)
                path.bypass(a, b)
            else:
                missing_paths.append(lineage)
        
        # Warn about missing paths
        if missing_paths:
            print(f"Warning: Paths not found in graph: {missing_paths}")

    def add_tdup(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds a tandem duplication in the specified paths, between positions a and b.
        Creates a loop in each path to represent the duplication.
        
        Parameters:
        - a (int): Starting position of the duplication
        - b (int): Ending position of the duplication
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        # Validate input parameters
        if a < 0 or b < 0: raise MSerror("Duplication positions must be positive")
        if a > b: raise MSerror("End must be > start")
        
        # Handle single lineage or collection
        lineages_to_process = [affected_lineages] if isinstance(affected_lineages, (int, str)) else affected_lineages
        missing_paths = []
        
        # Apply tandem duplication to each specified path
        for lineage in lineages_to_process:
            if lineage in self.paths:
                path = self.paths[lineage]
                # Check path length before operation
                if b > path.node_count:
                    raise MSerror(f"Path {lineage} too short for duplication at positions {a}-{b}")
                # Apply loop operation (creates tandem duplication)
                path.loop(a, b)
            else:
                missing_paths.append(lineage)
        
        # Warn about missing paths
        if missing_paths:
            print(f"Warning: Paths not found in graph: {missing_paths}")

    def add_inv(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds an inversion in the specified paths, between positions a and b.
        Creates an inverted section in each path to represent the inversion.
        
        Parameters:
        - a (int): Starting position of the inversion
        - b (int): Ending position of the inversion
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        # Validate input parameters
        if a <= 0 or b <= 0: raise MSerror("Inversion positions must be positive")
        if a > b: raise MSerror("End must be >= start")
        
        # Handle single lineage or collection
        lineages_to_process = [affected_lineages] if isinstance(affected_lineages, (int, str)) else affected_lineages
        missing_paths = []
        
        # Apply inversion to each specified path
        for lineage in lineages_to_process:
            if lineage in self.paths:
                path = self.paths[lineage]
                # Check path length before operation
                if b + 1 > path.node_count:
                    raise MSerror(f"Path {lineage} too short for inversion at positions {a}-{b}")
                # Apply invert operation (creates inversion)
                path.invert(a, b)
            else:
                missing_paths.append(lineage)
        
        # Warn about missing paths
        if missing_paths:
            print(f"Warning: Paths not found in graph: {missing_paths}")

if __name__ == "__main__":
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence("ABCDEFGHIJK",[1,2,3])
    graphA.details
    graphA.add_inv(2, 4, {3})
    graphA.details
    graphA.add_inv(1, 6, {3})
    graphA.details
    
    """
    print("start")
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence("ABCDEFG",[1,2,3])
    #graphA.details
    graphA.add_tdup(4, 4, {1, 2})
    #graphA.details

    nodeA = Node("|Je|", 0)
    nodeB = Node("|vais|", 1)
    nodeC = Node("|au|", 2)
    nodeD = Node("|concert|", 3)
    nodeE = Node("|ce|", 4)
    nodeF = Node("|soir|", 5)

    path = Path(1, [Edge(nodeA, True, nodeB, False),
                    Edge(nodeB, True, nodeC, True),
                    Edge(nodeC, False, nodeD, False),
                    Edge(nodeD, True, nodeE, False),
                    Edge(nodeE, True, nodeF, False)])
    print(path)
    print(repr(path))
    path.invert(start = 1, end = 3)
    print(path)
    print(repr(path))
    """

'''
class Graph:
    def add_deletion(self, start_idx: int, end_idx: int) -> None:
        deletion = DEL(A=self.nodes[start_idx], D=self.nodes[end_idx])
        deletion.compute_alt_seq()

    def add_insertion(self, idx: int, length) -> None:
        insertion = INS(A=self.nodes[idx], D=self.nodes[idx + 1], length=length)
        inserted_nodes = insertion.compute_alt_seq(self._node_id_generator)
        self.nodes = self.nodes + inserted_nodes 
    
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

