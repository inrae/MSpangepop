"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

This script constructs a variation graph from genomic sequences and mutation data,
representing structural variations across multiple lineages in a unified graph structure.
The resulting graph captures SNPs, insertions, deletions, inversions, and duplications
while maintaining the relationships between different evolutionary lineages.

Workflow:
    The main() function orchestrates the following pipeline:
    
    1. INITIALIZATION PHASE:
       - Create tracking objects for mutations (MutationRecap), graph optimization (LintVisualizer),
         and variant visualization (VariantSizeVisualizer)
       - Read input FASTA sequences (one per locus) and augmented traversal JSON
       - Calculate reference genome length from input sequences
    
    2. PARALLEL GRAPH CONSTRUCTION:
       - Extract lineage information from traversal data
       - Initialize individual graphs for each locus in parallel
       - Each graph represents ancestral sequence with paths for each lineage
       - Node creation is deferred (no IDs assigned yet)
    
    3. PARALLEL MUTATION APPLICATION:
       - Apply mutations from traversal data to corresponding graphs
       - Mutations include: SNPs, insertions (INS), deletions (DEL), 
         inversions (INV), duplications (DUP), replacements (REPL)
       - Track success/failure of each mutation application
       - Validate mutation positions against current path lengths
    
    4. GRAPH MERGING AND OPTIMIZATION:
       - Create GraphEnsemble container and concatenate all subgraphs
       - Calculate final sequence lengths for each lineage
       - Perform graph linting to remove orphan nodes (nodes not in any path)
       - Exclude ancestral path from final output if specified
    
    5. OUTPUT GENERATION:
       - Assign node IDs just-in-time during GFA serialization
       - Save graph in GFA v1.1 format with parallel I/O optimization
       - Export individual lineage sequences as FASTA files
       - Generate comprehensive visualization plots
       - Write detailed mutation recap file with statistics
    
REQUIRED INPUTS:
    --splited_fasta       : Path to split FASTA file containing sequence segments
                           (one sequence per tree in the ARG)
    --augmented_traversal : Path to JSON file containing the ORDERED mutation list and lineages
                           information from ARG traversal and mutation augmentation
"""

import itertools
import argparse
import psutil # type: ignore
from bitarray import bitarray # type: ignore
import threading
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning
from graph_utils import (
    mutate_base, generate_sequence, gather_lineages, 
    MutationRecap, VariantSizeVisualizer, LintVisualizer,
    GraphInitResult, MutationResult, ProgressTracker
)

# You can choose a matrix here for the SNP and insertion sequences
from matrix import random_matrix as snp_matrix, simple_at_bias_matrix as insertion_matrix

class Node:
    """Represents a node in the graph."""
    def __init__(self, dna_sequence, node_id=None): # Node ids allocation will append while saving
        self.id: int = node_id
        if isinstance(dna_sequence, str):
            self.dna_bases: bitarray = bitarray()
            self.dna_bases.frombytes(dna_sequence.encode("utf-8"))  # Encode string as bitarray
        else:
            self.dna_bases: bitarray = dna_sequence  # Store bitarray directly
            
    @property
    def __decode(self) -> str:
        """Decodes the bitarray back into a DNA sequence string."""
        return self.dna_bases.tobytes().decode("utf-8")

    def __repr__(self) -> str:
        """Returns only the raw DNA sequence as a string."""
        return self.__decode # Using str(Node) will automatically decode the node

    @property  
    def reverse_complement(self) -> str:
        """Returns the reverse complement of the DNA sequence."""
        sequence = self.__decode
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement_map.get(base, base) for base in reversed(sequence))


class Edge:
    """
    Represents a directed edge between two nodes, keeping track of the side of each node involved.
    
    Examples : 
        - Edge(A, +, B, -) is a direct link between A and B 
        - Edge(A, +, C, +) is a direct link between A and C but C will be read in reverse

    ODGI and VG call this class a "hook"
    """

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
        # Use object id if node.id is not set
        node1_id = self.node1.id if self.node1.id is not None else f"<{id(self.node1)}>"
        node2_id = self.node2.id if self.node2.id is not None else f"<{id(self.node2)}>"
        
        side1 = "+" if self.node1_side else "-"
        side2 = "+" if self.node2_side else "-"
        return f"({node1_id}{side1} -> {node2_id}{side2})"
    
class Path:
    """
    Represents a path in the graph, consisting of a ordered list of edges.

    Examples : 
        - Path([Edge(AB, +, CD, -), Edge(CD, +, EF, -)]) will produce the sequence AB CD EF
        - Path([Edge(AB, +, CD, +), Edge(CD, -, EF, -)]) will produce the sequence AB DC EF

    Each path, is in our graph a lineage.

    This class overloads the __getitem__(self, i) function so it yealds the node at the i position and NOT the Edge at i. 

    With edges [AB -> BC -> CD], we have:
        - 3 edges (node_count = 3)
        - 4 nodes (A, B, C, D)
        - Valid indices: 0, 1, 2, 3

    Using this we can do operation in the path sutch as : 
        - bypass
        - loop
        - invert
        - swap
        - paste
    between two nodes ids in the path. 
    """

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
            path_repr += f"{edge.node2 if not edge.node2_side else edge.node2.reverse_complement}"

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
            
            # Determine proper orientations for the loop-back edge
            # We need to exit from node[end] the opposite way we enter it
            # and enter node[start] the opposite way we exit it
            if end > 0 and end <= self.node_count:
                # Get how we enter the end node (from the edge coming into it)
                if end > 0 and end - 1 < len(self.path_edges):
                    end_entry_side = self.path_edges[end-1].node2_side
                    # Exit the opposite way
                    end_exit_side = not end_entry_side
                else:
                    end_exit_side = True  # Default if we can't determine
                
                # Get how we exit the start node
                if len(self.path_edges) > 0:
                    start_exit_side = self.path_edges[0].node1_side
                    # Enter the opposite way
                    start_entry_side = not start_exit_side
                else:
                    start_entry_side = False  # Default if we can't determine
            else:
                # Use defaults if indices are at boundaries
                end_exit_side = True
                start_entry_side = False
            
            # Create the "loop back" edge with proper orientations
            loop_back_edge = Edge(
                self[end],          # From end node
                end_exit_side,      # Exit orientation based on how we entered
                self[start],        # To start node (node 0)
                start_entry_side    # Entry orientation based on how we exit
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
            # Normal case: determine orientations based on existing edges
            # Internal edges within the duplicated section
            edges_to_duplicate = self.path_edges[start:end]
            
            # Determine proper orientations for the loop-back edge
            # based on how we're traversing these nodes in the path
            if end > 0 and end - 1 < len(self.path_edges):
                # How we enter the end node
                end_entry_side = self.path_edges[end-1].node2_side
                # Exit the opposite way for the loop-back
                end_exit_side = not end_entry_side
            else:
                end_exit_side = True  # Default
            
            if start > 0 and start < len(self.path_edges):
                # How we exit the start node  
                start_exit_side = self.path_edges[start].node1_side
                # Enter the opposite way for the loop-back
                start_entry_side = not start_exit_side
            else:
                start_entry_side = False  # Default
            
            # Create the "loop back" edge with proper orientations
            loop_back_edge = Edge(
                self[end],          # From end node
                end_exit_side,      # Proper exit orientation
                self[start],        # To start node
                start_entry_side    # Proper entry orientation
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
    """
    Represents a directed graph
    
    This graph holds : 
        - A set of unique nodes, that have no formal id and how a sequence (1 base for now)
        - Paths that reference an ordered list of Edges. Paths can share nodes but Edges are deepcopies.
        - The first and last node of the graph that cant be modified and are shared by all Paths. 
    """

    def __init__(self): 
        self.nodes = set()
        self.paths = {}
        self.start_node = None
        self.end_node = None
    
    def add_new_node(self, dna_sequence: str) -> Node:
        """Creates and adds a new node to the graph, return it"""
        node = Node(dna_sequence)
        self.nodes.add(node)
        return node

    def add_node(self, node: Node) -> None:
        """Adds a new node to the graph."""
        self.nodes.add(node)

    def __iadd__(self, other: "Graph") -> "Graph":
        """Merges another graph into this one. This combine the nodes and links similar paths"""

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
    
    def lint(self, ignore_ancestral=False, visualizer=None) -> None:
        """
        Removes orphan nodes from the graph.
        
        Orphan nodes are nodes that exist in self.nodes but are not referenced
        by any Edge in any Path. This can happen after graph modifications like
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

        self.nodes &= used_nodes # Use set operation for efficient linting 

    def build_from_sequence(self, nucleotide_sequence: str, lineages: set) -> None:
        """Constructs a graph and creates the associated paths from a nucleotide sequence."""
        if not nucleotide_sequence:
            raise MSerror("Cannot build a graph from an empty sequence.")
        
        ancestral_path = Path("Ancestral")
        
        # Create first node without ID
        self.start_node = Node(nucleotide_sequence[0])
        previous_node = self.start_node
        self.add_node(self.start_node)
        
        for base in nucleotide_sequence[1:]:
            current_node = self.add_new_node(base)
            ancestral_path += Edge(previous_node, True, current_node, False)
            previous_node = current_node
        
        self.end_node = previous_node
        self.paths[ancestral_path.lineage] = ancestral_path
        
        # Create deep copies of edges for each lineage
        for i in lineages:
            copied_path = Path(i)
            for edge in ancestral_path.path_edges:
                new_edge = Edge(
                    edge.node1,
                    edge.node1_side,
                    edge.node2,
                    edge.node2_side
                )
                copied_path.add_edge(new_edge)
            self.paths[i] = copied_path

    def __repr__(self) -> str:
        return f"Graph(id={self.id}, Nodes={len(self.nodes)}, Paths={len(self.paths)})"

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
    """
    This warper class is used to hold multiple Graph class. 

    We do this to handle each subgraph (locus) in paralell. 
    """

    def __init__(self, name: str = None, graph_list: list[Graph] = None):
        self.graphs: list[Graph] = graph_list if graph_list else []
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
                concatenated_graph += graph # We call the __iadd__ method of Graph
            self.graphs = [concatenated_graph]

    def lint(self, ignore_ancestral = False, visualizer=None): 
        """
        Warper method to lint all graphs in self, remove all node not in the paths
        """
        for graph in self: 
            graph.lint(ignore_ancestral, visualizer=visualizer)

# ============================================================================
# PARALLELIZATION FUNCTIONS
# ============================================================================

def parallel_graph_initialization(sequences, tree_lineages, chromosome, max_workers=4):
    """Parallelize graph initialization"""
    
    MScompute(f"Starting graph initialization for chromosome {chromosome}")

    progress = ProgressTracker(len(sequences), f"Initializing chromosome {chromosome} subgraphs")
    
    def init_single_graph(args):
        idx, record, tree_index, lineages = args
        
        sequence = str(record.seq)
        graph = Graph()
        graph.build_from_sequence(sequence, lineages)
        
        progress.update()
        
        return GraphInitResult(
            index=idx,
            graph=graph,
            tree_index=tree_index,
            lineages=lineages,
            node_count=len(graph.nodes)
        )
    
    # Prepare work items
    work_items = [(i, record, tree_index, lineages) 
                  for i, (record, (tree_index, lineages)) in enumerate(zip(sequences, tree_lineages))]
    
    graphs = [None] * len(sequences)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(init_single_graph, args): args[0] 
                  for args in work_items}
        
        for future in as_completed(futures):
            try:
                result = future.result()
                graphs[result.index] = result.graph
            except Exception as e:
                idx = futures[future]
                raise MSerror(f"Failed to initialize graph {idx}: {e}")
    
    return graphs

def parallel_apply_mutations(graphs, traversal, chromosome, max_workers=4):
    """
    Parallelize mutation application with thread-safe tracking
    """
    MScompute(f"Starting parallel mutation application for chromosome {chromosome}")
    
    # Progress tracker
    progress = ProgressTracker(len(graphs), f"Applying mutations to chromosome {chromosome} subgraphs")
    
    # Worker function for applying mutations to a single graph
    def apply_mutations_to_single_graph(args):
        graph_idx, graph, tree_data = args
        
        # Local tracking for this graph's mutations
        local_mutations = []
        local_variants = []
        
        # Extract tree metadata
        tree_index = tree_data.get("tree_index", "unknown")
        tree_interval = tree_data.get("initial_tree_interval", [0, 0])
        tree_start = tree_interval[0]
        
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
                continue
            
            # Apply each mutation
            for mutation in mutations:
                mut_type = mutation.get("type")
                start = mutation.get("start")
                length = mutation.get("length")
                
                # Handle None mutations
                if mut_type is None:
                    local_mutations.append({
                        "tree_index": tree_index,
                        "node_id": node_id,
                        "mutation_type": "SKIPPED",
                        "position": start,
                        "length": length,
                        "affected_lineages": affected_lineages,
                        "success": False,
                        "error_msg": "Mutation skipped in previous step"
                    })
                    continue
                
                if start is None:
                    local_mutations.append({
                        "tree_index": tree_index,
                        "node_id": node_id,
                        "mutation_type": mut_type,
                        "position": start,
                        "length": length,
                        "affected_lineages": affected_lineages,
                        "success": False,
                        "error_msg": "Mutation has no start position"
                    })
                    local_variants.append((mut_type, start, length, False, affected_lineages))
                    continue
                
                # Convert to relative position
                relative_start = start - tree_start
                
                # Validate and apply mutation
                success = False
                error_msg = None
                
                # Validate against each affected lineage's current path length
                failed_lineages = []
                error_messages = []
                
                for lineage in affected_lineages:
                    # Get the path for this lineage
                    path = graph.paths.get(lineage)
                    if not path:
                        failed_lineages.append(lineage)
                        error_messages.append(f"Lineage {lineage} not found in graph")
                        continue
                    
                    # Get current path length
                    current_path_length = path.node_count + 1
                    
                    # Validate based on mutation type
                    if mut_type == "SNP":
                        if relative_start == 0:
                            failed_lineages.append(lineage)
                            error_messages.append("SNP at first position would break connectivity")
                        elif relative_start >= current_path_length - 1:
                            failed_lineages.append(lineage)
                            error_messages.append(
                                f"SNP at or beyond last position (pos {relative_start}, path length {current_path_length})"
                            )
                    
                    elif mut_type == "INS":
                        if relative_start == 0:
                            failed_lineages.append(lineage)
                            error_messages.append("INS at position 0 would break connectivity")
                        elif relative_start >= current_path_length:
                            failed_lineages.append(lineage)
                            error_messages.append(
                                f"INS position {relative_start} >= path length {current_path_length}"
                            )
                    
                    elif mut_type == "DEL":
                        if length is None:
                            failed_lineages.append(lineage)
                            error_messages.append("DEL has no length specified")
                        else:
                            if relative_start == 0:
                                failed_lineages.append(lineage)
                                error_messages.append("DEL at first position would break connectivity")
                            elif relative_start + length > current_path_length - 1:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"DEL would affect last position (end {relative_start + length}, path length {current_path_length})"
                                )
                            elif length >= current_path_length - 2:
                                failed_lineages.append(lineage)
                                error_messages.append(
                                    f"DEL would leave less than 3 nodes (deleting {length} from {current_path_length} nodes)"
                                )
                    
                    elif mut_type == "INV":
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
                    valid_lineages = affected_lineages - set(failed_lineages)
                    error_msg = f"Invalid for lineages {failed_lineages}: {'; '.join(error_messages)}"
                    
                    # Record the failure
                    local_mutations.append({
                        "tree_index": tree_index,
                        "node_id": node_id,
                        "mutation_type": mut_type,
                        "position": start,
                        "length": length,
                        "affected_lineages": set(failed_lineages),
                        "success": False,
                        "error_msg": error_msg
                    })
                    
                    if not valid_lineages:
                        local_variants.append((mut_type, start, length, False, affected_lineages))
                        continue
                    
                    affected_lineages = valid_lineages
                
                # Try to apply the mutation
                try:
                    # Apply the mutation based on type
                    if mut_type == "SNP":
                        graph.add_snp(relative_start, affected_lineages)
                    elif mut_type == "INS":
                        if length is None:
                            raise MSerror(f"Insertion has no length")
                        graph.add_ins(relative_start, length, affected_lineages)
                    elif mut_type == "DEL":
                        if length is None:
                            raise MSerror(f"Deletion has no length")
                        graph.add_del(relative_start, relative_start + length, affected_lineages)
                    elif mut_type == "INV":
                        if length is None:
                            raise MSerror(f"Inversion has no length")
                        graph.add_inv(relative_start, relative_start + length - 1, affected_lineages)
                    elif mut_type == "DUP":
                        if length is None:
                            raise MSerror(f"Duplication has no length")
                        graph.add_tdup(relative_start, relative_start + length - 1, affected_lineages)
                    else:
                        raise MSerror(f"Unknown mutation type: {mut_type}")
                    
                    success = True
                    
                except Exception as e:
                    error_msg = str(e)
                    success = False
                
                # Track the mutation
                local_mutations.append({
                    "tree_index": tree_index,
                    "node_id": node_id,
                    "mutation_type": mut_type,
                    "position": start,
                    "length": length,
                    "affected_lineages": affected_lineages,
                    "success": success,
                    "error_msg": error_msg
                })
                
                # Track the variant
                local_variants.append((mut_type, start, length, success, affected_lineages))
        
        # Update progress
        progress.update()
        
        return MutationResult(
            index=graph_idx,
            mutations_applied=local_mutations,
            variants_tracked=local_variants,
            success=True
        )
    
    # Prepare work items
    work_items = [(i, graph, tree_data) 
                  for i, (graph, tree_data) in enumerate(zip(graphs, traversal))]
    
    # Collect results
    all_mutations = []
    all_variants = []
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(apply_mutations_to_single_graph, item): item[0] 
                  for item in work_items}
        
        # Collect results in order
        results = [None] * len(graphs)
        for future in as_completed(futures):
            try:
                result = future.result()
                results[result.index] = result
            except Exception as e:
                idx = futures[future]
                raise MSerror(f"Failed to apply mutations to graph {idx}: {e}")
        
        # Merge results in order
        for result in results:
            if result:
                all_mutations.extend(result.mutations_applied)
                all_variants.extend(result.variants_tracked)
    
    MScompute(f"Mutation application complete. Total mutations processed: {len(all_mutations)}")
    return all_mutations, all_variants

def apply_mutations_to_graphs(graphs, traversal, recap: MutationRecap, visualizer: VariantSizeVisualizer, chromosome, threads=4):
    """
    Apply mutations from the augmented traversal to the corresponding graphs.
    Now uses parallel processing by default.
    """
    # Use parallel mutation application
    mutations, variants = parallel_apply_mutations(graphs, traversal, chromosome, max_workers=threads)
    
    # Merge tracking results into recap and visualizer
    for mutation_data in mutations:
        recap.add_mutation(**mutation_data)
    
    for variant_data in variants:
        visualizer.add_variant(*variant_data)

def main(splited_fasta: str, augmented_traversal: str, output_file: str, 
         sample: str, chromosome: str, fasta_folder: str, recap_file: str = None, 
         variant_plot_dir: str = None, threads = 1) -> None:
    """
    Main function for graph creation with recap and visualization.
    Uses parallel processing by default for improved performance.
    """
    
    # Initialize tracking objects
    recap = MutationRecap(sample, chromosome)
    lint_visualizer = LintVisualizer()
    
    # Read input data first to get reference length
    MScompute("Reading input files")
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
        
        # Parallel graph initialization
        graphs = parallel_graph_initialization(
            sequences, tree_lineages, chromosome, 
            max_workers=threads
        )
        
        # Parallel mutation application
        apply_mutations_to_graphs(graphs, traversal, recap, var_visualizer, chromosome, threads)
        
        # Continue with concatenation and saving (sequential)
        ensemble = GraphEnsemble(name=f"{sample}_chr_{chromosome}", graph_list=graphs)

        MScompute(f"Concatenating and linting graphs")
        ensemble.concatenate
        if ensemble.graphs:  # Should have exactly one graph after concatenation
            final_graph = ensemble.graphs[0]
            lineage_lengths = {}
            for lineage, path in final_graph.paths.items():
                # Calculate the actual sequence length of the path
                lineage_lengths[lineage] = path.node_count + 1  # edges + 1 = nodes
            var_visualizer.set_lineage_lengths(lineage_lengths)

        ensemble.lint(ignore_ancestral=True, visualizer=lint_visualizer)

        MSpangepopDataHandler.save_to_gfav1_1_hybrid(
            ensemble,
            file_path=output_file,
            ignore_ancestral=True,
            max_workers=threads,
            sample=sample,
            chromosome=chromosome
        )
        
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
            MScompute(f"Saving recap [1/2]")
        
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
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing (default: 4)")

    args = parser.parse_args()
    
    main(args.splited_fasta, args.augmented_traversal, args.output_file, 
         args.sample, args.chromosome, args.fasta_folder, args.recap_file, 
         args.variant_plot_dir, args.threads)