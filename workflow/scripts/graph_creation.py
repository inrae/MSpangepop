import itertools
import argparse
import random
from bitarray import bitarray # type: ignore
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning
from graph_utils import mutate_base, generate_sequence

# You can choose a matrix here for the SNP and insertion sequences
from matrix import random_matrix as snp_matrix, simple_at_bias_matrix as insertion_matrix

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
    
    @property
    def bandage_path(self):
        if not self.path_edges:
            print("")
            return 
        path_repr = f"{self.path_edges[0].node1.id},"

        for edge in self.path_edges:
            path_repr += f"{edge.node2.id},"
        print(path_repr)

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
        Now handles boundary cases: start=0 (first node) and end=self.node_count (last node).
        
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
        
        # Get the start and end nodes
        start_node = self[start]
        end_node = self[end]
        
        # Determine the sides for the new edge
        # Handle different boundary cases
        
        if start == 0 and end == self.node_count:
            # Special case: bypass entire path (first to last node)
            # Use default orientation: True -> False
            start_side = True
            end_side = False
            
        elif start == 0:
            # Starting from first node - no incoming edge exists
            # Use the outgoing side from the first edge, but for the end node use the incoming side
            start_side = True  # Default orientation for start
            end_side = self.path_edges[end-1].node2_side  # Side of end node from its incoming edge
            
        elif end == self.node_count:
            # Ending at last node - no outgoing edge exists  
            # Use the incoming side from the start node, but for the end use default
            start_side = self.path_edges[start].node1_side  # Side of start node from its outgoing edge
            end_side = False  # Default orientation for end
            
        else:
            # Normal case: both nodes have incoming and outgoing edges
            start_side = self.path_edges[start].node1_side    # Side of start node from its outgoing edge
            end_side = self.path_edges[end-1].node2_side      # Side of end node from its incoming edge
        
        # Create the new bypass edge
        new_edge = Edge(start_node, start_side, end_node, end_side)
        
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
        - loop(0, 2) on A-B-C-D creates A-B-C-A-B-C-D (nodes A, B, and C are traversed twice)
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
        
        # Handle special case: start == 0 (duplicating from first node)
        if start == 0:
            # Edges that traverse the nodes [0, end]
            edges_to_duplicate = self.path_edges[start:end]
            
            # Create the "loop back" edge from node[end] to node[0]
            loop_back_edge = Edge(
                self[end],      # From end node
                True,           # Standard exit orientation
                self[start],    # To start node (node 0)
                False           # Standard entry orientation
            )
            
            # Create copies of the original edges for the second traversal
            duplicated_edges = []
            for edge in edges_to_duplicate:
                new_edge = Edge(
                    edge.node1,      # Same start node
                    edge.node1_side, # Same start node side
                    edge.node2,      # Same end node  
                    edge.node2_side  # Same end node side
                )
                duplicated_edges.append(new_edge)
            
            # Insert: loop_back_edge + duplicated_traversal at position 'end'
            complete_loop = [loop_back_edge] + duplicated_edges
            self.path_edges[end:end] = complete_loop
            
        else:
            # FIXED: Normal case now creates proper connectivity
            # We need to duplicate the INTERNAL edges of the section [start, end]
            # and add a loop-back edge from end to start
            
            # Internal edges within the duplicated section
            edges_to_duplicate = self.path_edges[start:end]  # FIXED: was start-1:end
            
            # Create the "loop back" edge from node[end] to node[start]
            loop_back_edge = Edge(
                self[end],      # From end node
                True,           # Standard exit orientation  
                self[start],    # To start node
                False           # Standard entry orientation
            )
            
            # Create copies of the internal edges
            duplicated_edges = []
            for edge in edges_to_duplicate:
                new_edge = Edge(
                    edge.node1,      # Same start node
                    edge.node1_side, # Same start node side
                    edge.node2,      # Same end node  
                    edge.node2_side  # Same end node side
                )
                duplicated_edges.append(new_edge)
            
            # Insert: loop_back_edge + duplicated_internal_edges at position 'end'
            complete_loop = [loop_back_edge] + duplicated_edges
            self.path_edges[end:end] = complete_loop
        
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
        
        # Dont update node count (should remain the same for inversions)

    def swap(self, pos: int, new_node: Node) -> None: 
        """
        Swap a node in the path with another node, retaining the orientation.
        Handles first, middle, and last nodes appropriately.
        
        Parameters:
        - pos (int): Position of the node to swap (0-based indexing)
        - new_node (Node): New node to replace the existing one
        """
        # Validation
        if pos < 0 or pos > self.node_count:
            raise MSerror("Position out of bounds")
        if self.node_count == 0:
            raise MSerror("Cannot swap in empty path")
        
        # Handle first node (pos == 0)
        if pos == 0:
            after_edge = self.path_edges[pos]  # Edge going out of node[0]
            after_edge.node1 = new_node

        # Handle last node (pos == self.node_count)  
        elif pos == self.node_count:
            # Only modify the edge coming INTO the last node
            before_edge = self.path_edges[pos-1]  # Edge coming into the last node
            before_edge.node2 = new_node
        
        # Handle middle nodes (normal case)
        else:
            # Modify both the edge coming in and the edge going out
            before_edge = self.path_edges[pos-1]  # Edge coming into the node
            after_edge = self.path_edges[pos]     # Edge going out of the node
            
            before_edge.node2 = new_node  # Update target of incoming edge
            after_edge.node1 = new_node   # Update source of outgoing edge

        # Dont update node count (should remain the same)

    def paste(self, start: int, end: int, list_of_nodes: list[Node]):
        """
        Insert nodes between two adjacent positions without removing anything.
        Only works when start + 1 == end (adjacent nodes).
        
        Example: paste(1, 2, [E, F]) on A→B→C→D gives A→B→E→F→C→D
        
        Parameters:
        - start (int): Starting node index  
        - end (int): Ending node index (must be start + 1)
        - list_of_nodes (list[Node]): List of nodes to insert
        """
        # Validation
        if start + 1 != end:
            raise MSerror("paste() only works between adjacent nodes (start + 1 must equal end)")
        if start > self.node_count or end > self.node_count:
            raise MSerror("Start and end indices out of bounds")
        if not list_of_nodes:
            raise MSerror("Cannot paste empty list of nodes")
        
        # Get the original edge between the adjacent nodes
        original_edge = self.path_edges[start]  # Edge connecting node[start] → node[end]
        
        # Handle single node insertion
        if len(list_of_nodes) == 1:
            single_node = list_of_nodes[0]
            
            # Create two edges: node[start] → single_node → node[end]
            first_edge = Edge(original_edge.node1, original_edge.node1_side, single_node, False)
            second_edge = Edge(single_node, True, original_edge.node2, original_edge.node2_side)
            
            # Replace the original edge with the two new edges
            self.path_edges[start:start+1] = [first_edge, second_edge]
        
        # Handle multiple node insertion
        else:
            # Create internal edges connecting the nodes in the list
            internal_edges = []
            for i in range(len(list_of_nodes) - 1):
                current_node = list_of_nodes[i]
                next_node = list_of_nodes[i + 1]
                internal_edges.append(Edge(current_node, True, next_node, False))
            
            # Create connecting edges
            connect_in = Edge(original_edge.node1, original_edge.node1_side, list_of_nodes[0], False)
            connect_out = Edge(list_of_nodes[-1], True, original_edge.node2, original_edge.node2_side)
            
            # Build complete insertion: connect_in + internal_chain + connect_out
            complete_insertion = [connect_in] + internal_edges + [connect_out]
            
            # Replace the original edge with the complete insertion
            self.path_edges[start:start+1] = complete_insertion
        
        # Update node count
        self.node_count = len(self.path_edges)

    def cut_paste(self, start: int, end: int, list_of_nodes: list[Node]):
        """
        Replace nodes between start and end with new nodes.
        Removes everything between start and end, then inserts the new nodes.
        
        Example: cut_paste(1, 3, [E, F]) on A→B→C→D gives A→B→E→F→D
        
        Parameters:
        - start (int): Starting node index
        - end (int): Ending node index (must be > start + 1)  
        - list_of_nodes (list[Node]): List of nodes to replace with
        """
        # Validation
        if start <= 0 or end <= 0:
            raise MSerror("Start and end indices must be positive")
        if start + 1 >= end:
            raise MSerror("cut_paste() requires end > start + 1 (use paste() for adjacent nodes)")
        if start > self.node_count or end > self.node_count:
            raise MSerror("Start and end indices out of bounds")
        if not list_of_nodes:
            raise MSerror("Cannot cut_paste with empty list of nodes")
        if start == 0 or end >= self.node_count:
            raise MSerror("Cannot cut_paste at path boundaries")
        
        # Get the edges that connect to the section we're replacing
        before_edge = self.path_edges[start-1]  # Edge leading TO node[start]  
        after_edge = self.path_edges[end]       # Edge leading FROM node[end]
        
        # Handle single node replacement
        if len(list_of_nodes) == 1:
            single_node = list_of_nodes[0]
            
            # Create edges connecting through the replacement node
            first_edge = Edge(before_edge.node1, before_edge.node1_side, single_node, False)
            second_edge = Edge(single_node, True, after_edge.node2, after_edge.node2_side)
            
            # Replace the entire section with the two new edges
            self.path_edges[start-1:end+1] = [first_edge, second_edge]
        
        # Handle multiple node replacement
        else:
            # Create internal edges connecting the replacement nodes
            internal_edges = []
            for i in range(len(list_of_nodes) - 1):
                current_node = list_of_nodes[i]
                next_node = list_of_nodes[i + 1]
                internal_edges.append(Edge(current_node, True, next_node, False))
            
            # Create connecting edges
            connect_in = Edge(before_edge.node1, before_edge.node1_side, list_of_nodes[0], False)
            connect_out = Edge(list_of_nodes[-1], True, after_edge.node2, after_edge.node2_side)
            
            # Build complete replacement: connect_in + internal_chain + connect_out
            complete_replacement = [connect_in] + internal_edges + [connect_out]
            
            # Replace the entire section with the complete replacement
            self.path_edges[start-1:end+1] = complete_replacement
        
        # Update node count
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

    def lint(self, ignore_ancestral=False) -> None:
        """
        Removes orphan nodes from the graph.
        
        Orphan nodes are nodes that exist in self.nodes but are not referenced
        by any edge in any path. This can happen after graph modifications like
        deletions, bypasses, or other operations that leave unused nodes behind.
        
        Parameters:
        - ignore_ancestral (bool): If True, nodes used only in the ancestral path
                                will be considered orphans and removed
        """
        # Collect all nodes that are actually used in path edges
        used_nodes = set()
        
        for lineage, path in self.paths.items():
            # Skip ancestral path if ignore_ancestral is True
            if ignore_ancestral and lineage == "Ancestral":
                continue
                
            for edge in path.path_edges:
                # Direct access
                used_nodes.add(edge.node1)
                used_nodes.add(edge.node2)
        
        # Keep only used nodes (more efficient than removing orphans)
        self.nodes &= used_nodes
        
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
        
        # Create deep copies of edges for each lineage
        for i in lineages:
            copied_path = Path(i)
            
            # Create new Edge objects for each path (deep copy)
            for edge in ancestral_path.path_edges:
                new_edge = Edge(
                    edge.node1,      # Same nodes
                    edge.node1_side, # Same orientation
                    edge.node2,      # Same nodes
                    edge.node2_side  # Same orientation
                )
                copied_path.add_edge(new_edge)
            
            self.paths[i] = copied_path

    def __repr__(self) -> str:
        return f"Graph(id={self.id}, Nodes={len(self.nodes)}, Paths={len(self.paths)})"

    @property
    def details(self):
        """Show graph details"""
        print(self)
        print("Paths :")
        for path in self.paths.items() :
            print(path)

    def _apply_to_paths(self, affected_lineages, operation_func):
        """
        Helper function to apply an operation to multiple paths.
        Handles single/multiple lineages and missing path warnings.
        
        Parameters:
        - affected_lineages: Single lineage or collection of lineages
        - operation_func: Function to apply to each valid path
        
        Returns:
        - List of valid (lineage, path) tuples that were processed
        """
        # Handle single lineage or collection
        lineages_to_process = [affected_lineages] if isinstance(affected_lineages, (int, str)) else affected_lineages
        missing_paths = []
        valid_paths = []
        
        # Collect valid paths and track missing ones
        for lineage in lineages_to_process:
            if lineage in self.paths:
                valid_paths.append((lineage, self.paths[lineage]))
            else:
                missing_paths.append(lineage)
        
        # Apply operation to each valid path
        for lineage, path in valid_paths:
            operation_func(path)
        
        # Warn about missing paths
        if missing_paths:
            MSwarning(f"Warning: Paths not found in graph: {missing_paths}")
        
        return valid_paths

    def add_del(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds a deletion between positions a and b by creating a bypass in each path.
        
        Parameters:
        - a (int): Starting position of the deletion
        - b (int): Ending position of the deletion  
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        self._apply_to_paths(affected_lineages, lambda path: path.bypass(a, b))

    def add_tdup(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds a tandem duplication between positions a and b by creating a loop in each path.
        
        Parameters:
        - a (int): Starting position of the duplication
        - b (int): Ending position of the duplication
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        self._apply_to_paths(affected_lineages, lambda path: path.loop(a, b))

    def add_inv(self, a: int, b: int, affected_lineages) -> None:
        """
        Adds an inversion between positions a and b by inverting the section in each path.
        
        Parameters:
        - a (int): Starting position of the inversion
        - b (int): Ending position of the inversion
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        self._apply_to_paths(affected_lineages, lambda path: path.invert(a, b))

    def add_snp(self, a: int, affected_lineages) -> None:
        """
        Adds a SNP at position a. Creates one shared mutated node across all affected paths.
        
        Parameters:
        - a (int): Position where the SNP should be added
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        # Get valid paths
        valid_paths = self._apply_to_paths(affected_lineages, lambda path: None)  # No operation yet
        
        if not valid_paths:
            return  # No valid paths to process
        
        # Get the original base from the first valid path
        first_lineage, first_path = valid_paths[0]
        node_base = str(first_path[a])
        
        # Verify consistency across all paths (sanity check)
        for lineage, path in valid_paths:
            if str(path[a]) != node_base:
                raise MSerror(f"Paths have different bases at position {a}: expected '{node_base}', found '{str(path[a])}' in path {lineage}")
        
        # Create one shared mutated node
        mutated_base = mutate_base(node_base, snp_matrix)
        shared_mutated_node = self.add_new_node(mutated_base)
        
        # Apply the shared node to all valid paths
        for lineage, path in valid_paths:
            path.swap(a, shared_mutated_node)

    def add_ins(self, a: int, length: int, affected_lineages) -> None:
        """
        Adds an insertion of specified length at position a using HMM-generated sequence.
        
        Parameters:
        - a (int): Position where the insertion should be added
        - length (int): Length of sequence to insert
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        if length <= 0:
            raise MSerror("Insertion length must be positive")
        
        # Generate insertion sequence using HMM
        insertion_sequence = generate_sequence(length, insertion_matrix)
        
        # Create nodes for the insertion sequence
        insertion_nodes = []
        for base in insertion_sequence:
            insertion_nodes.append(self.add_new_node(base))
        
        # Apply insertion to each valid path
        self._apply_to_paths(affected_lineages, lambda path: path.paste(a, a+1, insertion_nodes))

    def add_replacement(self, a: int, b: int, length: int, affected_lineages) -> None:
        """
        Replaces sequence between positions a and b with HMM-generated sequence of specified length.
        
        Parameters:
        - a (int): Starting position of the replacement
        - b (int): Ending position of the replacement
        - length (int): Length of replacement sequence
        - affected_lineages: Single lineage or collection of lineages to modify
        """
        if length <= 0:
            raise MSerror("Replacement length must be positive")
        
        # Generate replacement sequence using HMM
        replacement_sequence = generate_sequence(length, insertion_matrix)
        
        # Create nodes for the replacement sequence
        replacement_nodes = []
        for base in replacement_sequence:
            replacement_nodes.append(self.add_new_node(base))
        
        # Apply replacement to each valid path
        self._apply_to_paths(affected_lineages, lambda path: path.cut_paste(a, b, replacement_nodes))

class GraphEnsemble:
    def __init__(self, name: str = None, graph_list: list[Graph] = None):
        """
        Manages multiple graphs and allows linking them together.
        """
        self.graphs: list[Graph] = graph_list if graph_list else []
        self._node_id_generator: itertools.count = itertools.count(1)
        self.name = name
      
    def __repr__(self) -> str:
        if not self.graphs:
            return "GraphEnsemble(empty)"
        total_nodes = sum(len(graph.nodes) for graph in self.graphs)
        lines = [f"GraphEnsemble({len(self.graphs)} graphs, {total_nodes} total nodes"]
        for i, graph in enumerate(self.graphs):
            lines.append(f"\t[{i}] {graph}")
        return "\n".join(lines)
    
    def __getitem__(self, i: int) -> Graph: 
        """Returns the Graph at position i in the Ensemble"""
        if i > len(self.graphs): raise MSerror(f"Graph index {i} out of range, {len(self.graphs)}")
        return self.graphs[i]

    def add_graph(self, graph: Graph) -> None:
        self.graphs.append(graph)
    
    @property
    def concatenate(self) -> None:
            """
            Concatenates all graphs in the ensemble in place.
            Replaces the graph list with a single merged graph.
            """
            if not self.graphs:
                raise MSerror("No graphs to concatenate.")
            
            if len(self.graphs) == 1:
                return  # Already concatenated
            
            # Start with the first graph
            concatenated_graph = self.graphs[0]
            
            # Concatenate each subsequent graph
            for graph in self.graphs[1:]:
                concatenated_graph += graph
            self.graphs = [concatenated_graph]

    def lint(self, ignore_ancestral = False): 
        """
        Warper method to lint all graphs in self, remove all node not in the paths
        """
        for graph in self: 
            graph.lint(ignore_ancestral)


    # Care, method below is NOT AT ALL optimised, we should use node.outgoing_edges instead
    # But since we dont update them yet, we cant. 
    def save_to_gfav1_1(self, file_path, ignore_ancestral=False):
        """
        Saves the given graph to a GFA V1.1 (Graphical Fragment Assembly) file.
        The GFA format consists of:
        - "H" lines defining headers.
        - "S" lines defining sequence segments (nodes).
        - "L" lines defining links (edges) between segments.
        
        Parameters:
        - file_path (str): The output GFA file path.
        - ignore_ancestral (bool): If True, don't include edges from ancestral path
        
        Returns:
        None (writes to a file).
        """
        with open(file_path, 'w') as f:
            # Step 1: Write header with graph name
            f.write(f"H\t{self.name}\n")
            
            # Step 2: Write all nodes as sequence segments
            # Format: S\t{node_id}\t{sequence}\n
            for graph in self:
                for node in graph.nodes:
                    f.write(f"S\t{node.id}\t{node}\n")
            
                # Step 3: Collect all unique edges from paths
                unique_edges = set()
                
                for lineage, path in graph.paths.items():
                    # Skip ancestral path if ignore_ancestral is True
                    if ignore_ancestral and lineage == "Ancestral":
                        continue
                    
                    # Add all edges from this path to the set
                    for edge in path.path_edges:
                        # Create a tuple that uniquely identifies this edge
                        # (node1_id, node1_side, node2_id, node2_side)
                        edge_signature = (
                            edge.node1.id,
                            edge.node1_side,
                            edge.node2.id,
                            edge.node2_side
                        )
                        unique_edges.add(edge_signature)
                
                # Step 4: Write all unique edges as links
                # Format: L\t{node1_id}\t{orientation1}\t{node2_id}\t{orientation2}\t0M\n
                for node1_id, node1_side, node2_id, node2_side in unique_edges:
                    # Convert sides to GFA orientation symbols
                    # True = +, False = - for node1_side
                    orientation1 = "+" if node1_side else "-"
                    
                    # False = +, True = - for node2_side (as specified)
                    orientation2 = "+" if not node2_side else "-"
                    
                    f.write(f"L\t{node1_id}\t{orientation1}\t{node2_id}\t{orientation2}\t0M\n")
                
                # Step 5, Write the paths as W-lines
                for lineage, path in graph.paths.items():
                    # Skip ancestral path if ignore_ancestral is True
                    if ignore_ancestral and lineage == "Ancestral":
                        continue
                    f.write(f"W\tlineage_{path.lineage}\t0\tlineage_{path.lineage}\t0\t0\t>{repr(path)}\n")
                    

if __name__ == "__main__":
    count = itertools.count(1)
    graphA = Graph(count)
    graphA.build_from_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",[1,2,3,4])
    graphA.add_del(6,12, affected_lineages={1})
    graphA.add_ins(15 ,12, affected_lineages={1})
    graphA.add_snp(3, affected_lineages={1})
    graphA.add_inv(30 ,35, affected_lineages={1})
    graphA.add_tdup(50 ,60, affected_lineages={1})

    chromosome = GraphEnsemble(name = "Test Graph", graph_list = [graphA])
    chromosome.lint(ignore_ancestral=True)
    chromosome.save_to_gfav1_1("./test.gfa", ignore_ancestral=True)


'''
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
