"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : This scripts create a variation graph from a json file full of mutations and lineages.

--splited_fasta ath to split FASTA file
--augmented_traversal Path to augmented traversal JSON
--output_file Path to output GFA file
--sample Sample name
--chromosome Chromosome identifier
--recap_file Path to save recap file
--variant_plot_dir Directory to save variant size plots
"""

import itertools
import argparse
import psutil # type: ignore
from bitarray import bitarray # type: ignore
import threading
from io import StringIO
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
import time

from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning
from graph_utils import mutate_base, generate_sequence,gather_lineages, MutationRecap, VariantSizeVisualizer, LintVisualizer

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

        self.false_side: set["Edge"] = set()  # We use this for merging
        self.true_side: set["Edge"] = set()

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

    def __eq__(self, other):
        "Check if two edges are the same"
        if not isinstance(other, Edge):
            return NotImplemented
        return (
            self.node1 == other.node1 and
            self.node1_side == other.node1_side and
            self.node2 == other.node2 and
            self.node2_side == other.node2_side
        )

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
        
        Note: node_count represents the number of EDGES in the path.
        The actual number of nodes = node_count + 1 (edges + 1)
        """
        self.lineage = lineage
        self.path_edges: list[Edge] = path_edges if path_edges else []
        self.node_count = len(self.path_edges)  # This is actually edge count!

    @property
    def num_nodes(self) -> int:
        """Returns the actual number of nodes in the path."""
        if not self.path_edges:
            return 0
        return self.node_count + 1  # edges + 1 = nodes

    def add_edge(self, edge: Edge) -> None:
        """Adds an edge to the path."""
        self.path_edges.append(edge)
        self.node_count += 1  # Increment edge count

    def __repr__(self) -> str:
        """Returns the gfa string of the path"""
        if not self.path_edges:
            return ""  # If there are no edges, return an empty string

        # Start with the first node
        path_repr = f"{self.path_edges[0].node1.id}+"

        # Then for each edge, append the node2 id (only if it's not the first node)
        for edge in self.path_edges:
            path_repr += f",{edge.node2.id}{'+' if not edge.node2_side else '-'}"

        return path_repr

    def __str__(self) -> str:
        """Returns a string of the path"""
        if not self.path_edges:
            return ""

        path_repr = f"{self.path_edges[0].node1}"

        for edge in self.path_edges:
            path_repr += f"{edge.node2 if not edge.node2_side else edge.node2.reversed}"

        return path_repr

    def __getitem__(self, i: int) -> Node:
        """Returns the node at position i in the path, supporting negative indices.
        
        Important: With edges [AB -> BC -> CD], we have:
        - 3 edges (node_count = 3)
        - 4 nodes (A, B, C, D)
        - Valid indices: 0, 1, 2, 3
        """
        if not self.path_edges:
            raise MSerror("Path is empty")

        if i < 0:
            i = self.node_count + 1 + i  # Convert negative index to positive

        if i < self.node_count:
            return self.path_edges[i].node1
        elif i == self.node_count:
            return self.path_edges[-1].node2
        else:
            raise MSerror(f"Node index {i} out of bounds for path with {self.node_count + 1} nodes")

    def __iadd__(self, other):
        """Allow concatenation of two paths (or an edge to a path)"""
        if isinstance(other, Edge):
            self.add_edge(other)  # If we try to add an edge
        elif isinstance(other, Path):
            for path_edge in other.path_edges:
                self.add_edge(path_edge)  # If we try to add a path
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
        if start < 0 or end < 0: 
            raise MSerror("Start and end indices must be non-negative")
        if start >= self.node_count or end >= self.node_count: 
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

    def lint(self, ignore_ancestral=False, visualizer=None) -> None:
        """
        Removes orphan nodes from the graph.
        
        Orphan nodes are nodes that exist in self.nodes but are not referenced
        by any edge in any path. This can happen after graph modifications like
        deletions, bypasses, or other operations that leave unused nodes behind.
        
        Parameters:
        - ignore_ancestral (bool): If True, nodes used only in the ancestral path
                                will be considered orphans and removed
        """
        used_nodes = set()

        for lineage, path in self.paths.items():
            if ignore_ancestral and lineage == "Ancestral":
                continue
            for edge in path.path_edges:
                used_nodes.add(edge.node1)
                used_nodes.add(edge.node2)

        if visualizer:
            all_nodes = set(self.nodes)
            removed_nodes = all_nodes - used_nodes
            visualizer.record(
                before=len(all_nodes),
                after=len(used_nodes),
                removed=removed_nodes
            )

        self.nodes &= used_nodes
    
    def __iter__(self):
        """Make GraphEnsemble iterable."""
        return iter(self.graphs)

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

    def lint(self, ignore_ancestral = False, visualizer=None): 
        """
        Warper method to lint all graphs in self, remove all node not in the paths
        """
        for graph in self: 
            graph.lint(ignore_ancestral, visualizer=visualizer)

    def save_to_gfav1_1_hybrid(self, file_path, ignore_ancestral=False, max_workers=None):
        """
        Hybrid approach with progress reporting: parallel data collection + sequential writing.
        Best balance of performance and simplicity with detailed progress updates.
        """
        
        if not self.graphs:
            raise MSerror("No graphs to save")
        
        graph = self.graphs[0]
        start_time = time.time()
        
        # Set default max_workers if not provided
        if max_workers is None:
            max_workers = 4
        
        # Count total items for progress reporting
        total_nodes = len(graph.nodes)
        total_paths = sum(1 for lineage in graph.paths.keys() 
                        if not (ignore_ancestral and lineage == "Ancestral"))
        print()
        MScompute(f"Starting GFA export: {total_nodes} nodes, {total_paths} paths to process")
        
        def process_nodes_parallel():
            """Process nodes with chunking if beneficial"""

            nodes_list = list(graph.nodes)

            chunk_size = max(1, len(nodes_list) // max_workers)
            
            def process_chunk(chunk_data):
                chunk_idx, chunk = chunk_data
                buffer = StringIO()
                for node in chunk:
                    buffer.write(f"S\t{node.id}\t{node}\n")
                return chunk_idx, buffer.getvalue()
            
            # Fixed: Create chunks properly with actual node objects
            chunks = []
            for i in range(0, len(nodes_list), chunk_size):
                chunk_nodes = nodes_list[i:i + chunk_size]
                chunk_idx = len(chunks)
                chunks.append((chunk_idx, chunk_nodes))
            
            completed_chunks = 0
            results = [''] * len(chunks)
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_chunk = {executor.submit(process_chunk, chunk_data): chunk_data 
                                for chunk_data in chunks}
                
                for future in as_completed(future_to_chunk):
                    chunk_idx, chunk_result = future.result()
                    results[chunk_idx] = chunk_result
                    completed_chunks += 1
            
            result = ''.join(results)
    
            return result
        
        def process_paths_and_edges():
            """Process paths and collect edges simultaneously"""
            
            paths_buffer = StringIO()
            unique_edges = set()
            processed_paths = 0
            
            for lineage, path in graph.paths.items():
                if ignore_ancestral and lineage == "Ancestral":
                    continue
                
                # Progress update every 100 paths or for large datasets
                if processed_paths % max(1, total_paths // 10) == 0 and processed_paths > 0:
                    progress = (processed_paths / total_paths) * 100
                
                # Process path
                paths_buffer.write(f"P\tlineage_{path.lineage}\t{repr(path)}\n")
                
                # Collect edges from this path
                path_edges_count = 0
                for edge in path.path_edges:
                    edge_signature = (
                        edge.node1.id, edge.node1_side,
                        edge.node2.id, edge.node2_side
                    )
                    unique_edges.add(edge_signature)
                    path_edges_count += 1
                
                processed_paths += 1
            
            MScompute("Converting edges to GFA format...")
            
            # Convert edges to strings
            edges_buffer = StringIO()
            processed_edges = 0
            
            for edge_data in unique_edges:

            
                node1_id, node1_side, node2_id, node2_side = edge_data
                orientation1 = "+" if node1_side else "-"
                orientation2 = "+" if not node2_side else "-"
                edges_buffer.write(f"L\t{node1_id}\t{orientation1}\t{node2_id}\t{orientation2}\t0M\n")
                processed_edges += 1
            

            return paths_buffer.getvalue(), edges_buffer.getvalue()
        
        # Execute both tasks in parallel
        MScompute("Starting parallel processing of nodes and paths...")
        
        with ThreadPoolExecutor(max_workers=2) as executor:
            nodes_future = executor.submit(process_nodes_parallel)
            paths_edges_future = executor.submit(process_paths_and_edges)
            
            nodes_str = nodes_future.result()
            paths_str, edges_str = paths_edges_future.result()
        
        # Calculate data sizes for reporting
        nodes_size = len(nodes_str)
        edges_size = len(edges_str)
        paths_size = len(paths_str)
        total_size = nodes_size + edges_size + paths_size
        
        MScompute(f"Data prepared: {total_size:,} characters ({nodes_size:,} nodes, {edges_size:,} edges, {paths_size:,} paths)")
        
        # Write to file
        MScompute(f"Writing GFA file to {file_path}...")
        write_start = time.time()
        
        with open(file_path, 'w') as f:
            print("\t\tWriting nodes...")
            f.write(nodes_str)
            
            print("\t\tWriting edges...")
            f.write(edges_str)
            
            print("\t\tWriting paths...")
            f.write(paths_str)
        
        write_time = time.time() - write_start
        total_time = time.time() - start_time
        
        # Final statistics
        file_size_mb = total_size / (1024 * 1024)
        write_speed_mb_s = file_size_mb / write_time if write_time > 0 else 0
        print()
        MScompute(f"GFA export completed successfully!")
        print(f"\t\tSize: {file_size_mb:.2f} MB")
        print(f"\t\tWrite speed: {write_speed_mb_s:.2f} MB/s")
        print()

def apply_mutations_to_graphs(graphs, traversal, recap: MutationRecap, visualizer: VariantSizeVisualizer, chromosome):
    """
    Apply mutations from the augmented traversal to the corresponding graphs.
    
    Simplified version that only validates against actual path lengths, not static tree boundaries.
    """
    i = 0 
    # Process each tree and its corresponding graph
    for tree_idx, (tree_data, graph) in enumerate(zip(traversal, graphs)):
        # Extract tree metadata (for reporting only)
        tree_index = tree_data.get("tree_index", "unknown")
        tree_interval = tree_data.get("initial_tree_interval", [0, 0])
        tree_start = tree_interval[0]
        
        if i == 0 or i % 10 == 0 : 
            MScompute(f"Applying mutations to chromosome {chromosome} subgraphs -> {(tree_index*100)/len(graphs):.0f} %")
        
        # Process each node in the tree
        for node_data in tree_data.get("nodes", []):
            node_id = node_data.get("node")
            mutations = node_data.get("mutations", [])
            affected_nodes = set(node_data.get("affected_nodes", []))
            
            if not mutations:
                continue
            
            lineages = set(tree_data.get("lineages", []))
            affected_lineages = affected_nodes.intersection(lineages)
            
            if not affected_lineages:
                MSwarning(f"Node {node_id} has mutations but no affected lineages")
                continue
            
            # Apply each mutation for this node
            for mutation in mutations:
                mut_type = mutation.get("type")
                start = mutation.get("start")
                length = mutation.get("length")
                
                # HANDLE NONE MUTATIONS - Skip mutations that couldn't be placed
                if mut_type is None:
                    # This mutation was skipped in draw_variants.py due to interval constraints
                    error_msg = "Mutation skipped in previous step (locus got too small)"
                    recap.add_mutation(tree_index, node_id, "SKIPPED", start, length, 
                                     affected_lineages, False, error_msg)
                    continue
                
                # Check if start position is None
                if start is None:
                    error_msg = "Mutation has no start position"
                    recap.add_mutation(tree_index, node_id, mut_type, start, length, 
                                     affected_lineages, False, error_msg)
                    visualizer.add_variant(mut_type, start, length, False, affected_lineages)
                    continue
                
                # Convert to relative position within the tree
                relative_start = start - tree_start
                
                # Validate against each affected lineage's CURRENT path length
                # Path lengths can vary between lineages due to previous mutations
                failed_lineages = []
                error_messages = []
                
                for lineage in affected_lineages:
                    # Get the path for this lineage
                    path = graph.paths.get(lineage)
                    if not path:
                        failed_lineages.append(lineage)
                        error_messages.append(f"Lineage {lineage} not found in graph")
                        continue
                    
                    # Get current path length (number of nodes)
                    # CRITICAL FIX: node_count is edges, actual nodes = edges + 1
                    current_path_length = path.node_count + 1
                    
                    # Validate based on mutation type - only check against current path
                    # No need to check static tree boundaries!
                    
                    if mut_type == "SNP":
                        # SNPs modify a single existing position [0, path_length-1]
                        # Cannot be at first or last position for connectivity
                        if relative_start == 0:
                            failed_lineages.append(lineage)
                            error_messages.append("SNP at first position would break connectivity")
                        elif relative_start >= current_path_length - 1:
                            failed_lineages.append(lineage)
                            error_messages.append(
                                f"SNP at or beyond last position (pos {relative_start}, path length {current_path_length})"
                            )
                    
                    elif mut_type == "INS":
                        # Insertions add sequence between two adjacent positions
                        # Cannot insert at position 0 (before first node)
                        # Can insert up to position path_length-1 (after last node)
                        if relative_start == 0:
                            failed_lineages.append(lineage)
                            error_messages.append("INS at position 0 would break connectivity")
                        elif relative_start >= current_path_length:
                            failed_lineages.append(lineage)
                            error_messages.append(
                                f"INS position {relative_start} >= path length {current_path_length}"
                            )
                    
                    elif mut_type == "DEL":
                        # Deletion removes nodes from [start, start+length) 
                        # Cannot delete first or last node
                        if length is None:
                            failed_lineages.append(lineage)
                            error_messages.append("DEL has no length specified")
                        else:
                            # Check boundaries against current path
                            if relative_start == 0:
                                failed_lineages.append(lineage)
                                error_messages.append("DEL at first position would break connectivity")
                            elif relative_start + length > current_path_length - 1:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"DEL would affect last position (end {relative_start + length}, path length {current_path_length})"
                                )
                            # Also ensure we don't delete too much
                            elif length >= current_path_length - 2:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"DEL would leave less than 3 nodes (deleting {length} from {current_path_length} nodes)"
                                )
                    
                    elif mut_type == "INV":
                        # Inversion affects nodes [start, start+length-1] inclusive
                        # Cannot invert first or last node
                        if length is None:
                            failed_lineages.append(lineage)
                            error_messages.append("INV has no length specified")
                        else:
                            actual_end = relative_start + length - 1
                            
                            if relative_start == 0:
                                failed_lineages.append(lineage)
                                error_messages.append("INV at first position would break connectivity")
                            elif actual_end >= current_path_length - 1:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"INV would affect last position (end {actual_end}, path length {current_path_length})"
                                )
                    
                    elif mut_type == "DUP":
                        # Duplication affects nodes [start, start+length-1] inclusive
                        # Cannot duplicate first or last node
                        if length is None:
                            failed_lineages.append(lineage)
                            error_messages.append("DUP has no length specified")
                        else:
                            actual_end = relative_start + length - 1
                            
                            if relative_start == 0:
                                failed_lineages.append(lineage)
                                error_messages.append("DUP at first position would break connectivity")
                            elif actual_end >= current_path_length - 1:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"DUP would affect last position (end {actual_end}, path length {current_path_length})"
                                )
                
                # Handle validation failures
                if failed_lineages:
                    # Some lineages failed validation
                    valid_lineages = affected_lineages - set(failed_lineages)
                    error_msg = f"Invalid for lineages {failed_lineages}: {'; '.join(error_messages)}"
                    
                    # Record the failure for failed lineages
                    recap.add_mutation(tree_index, node_id, mut_type, start, length, 
                                     set(failed_lineages), False, error_msg)
                    
                    # If no valid lineages remain, skip this mutation entirely
                    if not valid_lineages:
                        visualizer.add_variant(mut_type, start, length, False, affected_lineages)
                        continue
                    
                    # Otherwise, continue with valid lineages only
                    affected_lineages = valid_lineages
                
                # Try to apply the mutation to valid lineages
                success = False
                try:
                    # Apply the appropriate mutation type
                    # Note: Graph methods use 0-based indexing
                    
                    if mut_type == "SNP":
                        # Single nucleotide polymorphism at one position
                        graph.add_snp(relative_start, affected_lineages)
                    
                    elif mut_type == "INS":
                        # Insertion of new sequence at a position
                        if length is None:
                            raise MSerror(f"Insertion at position {start} has no length specified")
                        graph.add_ins(relative_start, length, affected_lineages)
                    
                    elif mut_type == "DEL":
                        # Deletion from start to start+length (exclusive)
                        if length is None:
                            raise MSerror(f"Deletion at position {start} has no length specified")
                        graph.add_del(relative_start, relative_start + length, affected_lineages)
                    
                    elif mut_type == "INV":
                        # Inversion from start to start+length-1 (inclusive)
                        if length is None:
                            raise MSerror(f"Inversion at position {start} has no length specified")
                        graph.add_inv(relative_start, relative_start + length - 1, affected_lineages)
                    
                    elif mut_type == "DUP":
                        # Tandem duplication from start to start+length-1 (inclusive)
                        if length is None:
                            raise MSerror(f"Duplication at position {start} has no length specified")
                        graph.add_tdup(relative_start, relative_start + length - 1, affected_lineages)
                    
                    else:
                        # Unknown mutation type
                        raise MSerror(f"Unknown mutation type: {mut_type}")
                    
                    # Mutation applied successfully
                    success = True
                    
                    # Record success
                    recap.add_mutation(tree_index, node_id, mut_type, start, length, 
                                     affected_lineages, True)
                    
                except Exception as e:
                    # Mutation application failed
                    error_msg = str(e)
                    
                    # Record failure
                    recap.add_mutation(tree_index, node_id, mut_type, start, length, 
                                     affected_lineages, False, error_msg)
                
                # Track variant size for visualization (both success and failure)
                visualizer.add_variant(mut_type, start, length, success, affected_lineages)
        
        i += 1


def main(splited_fasta: str, augmented_traversal: str, output_file: str, 
         sample: str, chromosome: str, fasta_folder: str, recap_file: str = None, 
         variant_plot_dir: str = None, threads = 1) -> None:
    """Main function for graph creation with recap and visualization."""
    
    # Initialize tracking objects
    recap = MutationRecap(sample, chromosome)
    lint_visualizer = LintVisualizer()
    MScompute("Starting to generate the variation graph")
    
    # Read input data first to get reference length
    sequences = MSpangepopDataHandler.read_fasta(splited_fasta)
    traversal = MSpangepopDataHandler.read_json(augmented_traversal)
    
    # Calculate reference length (sum of all original sequence lengths)
    reference_length = sum(len(str(record.seq)) for record in sequences)
    
    # Initialize visualizer with reference length
    var_visualizer = VariantSizeVisualizer(sample, chromosome, reference_length)
    
    try:
        tree_lineages = gather_lineages(traversal)
        
        if len(sequences) != len(tree_lineages):
            raise MSerror(f"Mismatch: {len(sequences)} sequences but {len(tree_lineages)} trees")
        
        node_id_generator = itertools.count(1)
        
        graphs = []
        y = 0
        for i, (record, (tree_index, lineages)) in enumerate(zip(sequences, tree_lineages)):

            if y == 0 or y % 50 == 0 :
                memory = psutil.virtual_memory()
                available_gb = memory.available / (1024**3)
                MScompute(f"Initialization of chromosome {chromosome} subgraphs -> {(tree_index*100)/len(sequences):.0f}% | {available_gb:.1f} GB available")
                            
            sequence = str(record.seq)
            new_graph = Graph(node_id_generator)
            new_graph.build_from_sequence(sequence, lineages)
            graphs.append(new_graph)
            y += 1
        
        MScompute(f"Initialized {len(graphs)} graphs for chromosome {chromosome}")
        
        # Apply mutations with tracking

        MScompute(f"Starting to integrate mutation to all graphs")
        apply_mutations_to_graphs(graphs, traversal, recap, var_visualizer, chromosome)
        ensemble = GraphEnsemble(name=f"{sample}_chr_{chromosome}", graph_list=graphs)

        MScompute(f"Concatenating all graphs")
        ensemble.concatenate
        if ensemble.graphs:  # Should have exactly one graph after concatenation
            final_graph = ensemble.graphs[0]
            lineage_lengths = {}
            for lineage, path in final_graph.paths.items():
                # Calculate the actual sequence length of the path
                lineage_lengths[lineage] = path.node_count + 1  # edges + 1 = nodes
            var_visualizer.set_lineage_lengths(lineage_lengths)

        MScompute(f"Removing nodes orphaned by the mutations on {chromosome}")
        ensemble.lint(ignore_ancestral=True, visualizer=lint_visualizer)

        ensemble.save_to_gfav1_1_hybrid(output_file, ignore_ancestral=True, max_workers=threads)
        
        MSpangepopDataHandler.write_fasta_multithreaded(
            graph=final_graph,
            sample=sample,
            chromosome=chromosome,
            fasta_folder=fasta_folder,
            threads=threads
        )

    except Exception as e:
        recap.summary["fatal_error"] = str(e)
        raise
    
    finally:
        # Save recap
        if recap_file:
            recap.save_recap(recap_file)
            MScompute(f"Recap saved to {recap_file}")
        
        if variant_plot_dir:
            
            os.makedirs(variant_plot_dir, exist_ok=True)
            
            dist_plot_path = os.path.join(variant_plot_dir, 
                                        f"{sample}_chr{chromosome}_graph_variant_sizes.png")
            var_visualizer.save_size_distribution_plot(dist_plot_path)
            
            density_plot_path = os.path.join(variant_plot_dir,
                                            f"{sample}_chr{chromosome}_graph_variant_density.png")
            var_visualizer.save_variant_density_plot(density_plot_path)
            
            shared_plot_path = os.path.join(variant_plot_dir,
                                        f"{sample}_chr{chromosome}_graph_shared_variants.png")
            var_visualizer.save_shared_variants_heatmap(shared_plot_path)
            
            proportions_plot_path = os.path.join(variant_plot_dir,
                                            f"{sample}_chr{chromosome}_graph_variant_proportions.png")
            var_visualizer.save_variant_type_proportions_plot(proportions_plot_path)
            
            cumulative_plot_path = os.path.join(variant_plot_dir,
                                            f"{sample}_chr{chromosome}_graph_cumulative_variants.png")
            var_visualizer.save_cumulative_variants_plot(cumulative_plot_path)
            
            lengths_plot_path = os.path.join(variant_plot_dir,
                                        f"{sample}_chr{chromosome}_graph_lineage_lengths.png")
            var_visualizer.save_lineage_lengths_plot(lengths_plot_path)

            lint_visualizer.write_txt_report(recap_file)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create variation graph from mutations")
    parser.add_argument("--splited_fasta", required=True, help="Path to split FASTA file")
    parser.add_argument("--augmented_traversal", required=True, help="Path to augmented traversal JSON")
    parser.add_argument("--output_file", required=True, help="Path to output GFA file")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--chromosome", required=True, help="Chromosome identifier")
    parser.add_argument("--recap_file", help="Path to save recap file")
    parser.add_argument("--variant_plot_dir", help="Directory to save variant size plots")
    parser.add_argument("--fasta_folder", help="Directory to save all fasta")
    parser.add_argument("--threads", type=int)

    args = parser.parse_args()
    main(args.splited_fasta, args.augmented_traversal, args.output_file, 
         args.sample, args.chromosome,args.fasta_folder, args.recap_file, args.variant_plot_dir, args.threads)  