"""
================================================================================
                    NODE MERGING OPTIMIZATION ALGORITHM
================================================================================

Author: Lucien Piat (INRAe, PangenOak Project)
Based on insights from: Sigfried Dubois

Description:
    This algorithm performs graph simplification by merging nodes that have 
    identical connectivity patterns. Two nodes can be merged when:
    - All edges from node1's true side go to the same target (node2)
    - All edges to that side of node2 come exclusively from node1's true side
    
    Two merge types are supported:
    1. True→False merge: Prepends node1's sequence to node2
    2. True→True merge: Appends node1's reversed sequence to node2

Algorithm Complexity:
    Original Implementation: O(M × N × L × E)
    - M: Number of merges performed
    - N: Number of nodes in graph  
    - L: Number of lineages (paths)
    - E: Average edges per path
    - Worst case O(N² × L × E) when M = N
    
    Optimized Implementation: O(I × (N + L × E))
    - I: Number of iterations (typically log N in practice)
    - Builds adjacency data structures once per iteration: O(L × E)
    - Finds all candidates in single pass: O(N)
    - Processes non-conflicting merges: O(M × L × E) where M ≤ N
    
    Key optimizations:
    - Batch candidate detection instead of one-by-one
    - Process multiple merges per iteration
    - Efficient set-based conflict detection
    - No redundant graph traversals

Example:
    Initial: A→B→C and A→B→D with B having only incoming edges from A
    Result: A→BC and A→BD (B merged into C and D)

================================================================================
"""


from collections import defaultdict
from io_handler import MSerror, MSsuccess, MSwarning, MScompute
from collections import defaultdict
from bitarray import bitarray # type: ignore

def merge_node_True_False(node1, node2, graph):
    """
    Merge node1 INTO node2 by prepending node1's sequence to node2.
    After this, node1 will be removed and node2 will contain both sequences.
    """

    # Step 1: Prepend node1's bases to node2's sequence
    node2.dna_bases = node1.dna_bases + node2.dna_bases
    
    # Step 2: Update all paths
    for lineage, path in graph.paths.items():
        edges_to_remove = []
        
        for i, edge in enumerate(path.path_edges):
            # Remove the connecting edge between node1 and node2
            if edge.node1 == node1 and edge.node2 == node2:
                edges_to_remove.append(i)
                
            # Redirect ALL edges involving node1 to node2
            else:
                if edge.node1 == node1:
                    edge.node1 = node2
                if edge.node2 == node1:
                    edge.node2 = node2
        
        # Remove connecting edges
        for i in reversed(edges_to_remove):
            del path.path_edges[i]
        
        path.node_count = len(path.path_edges)

def merge_node_True_True(node1, node2, graph):
    """
    Merge node1 INTO node2 by appending node1's reverse sequence to node2.
    After this, node1 will be removed and node2 will contain both sequences.
    """
    # Step 1: Get the reverse sequence of node1 and append it to node2
    # Decode to string, reverse, then encode back to preserve UTF-8 encoding
    node1_sequence = node1.dna_bases.tobytes().decode("utf-8")
    node1_reversed = node1_sequence[::-1]  # Reverse at character level
    
    # Convert back to bitarray
    node1_reversed_bits = bitarray()
    node1_reversed_bits.frombytes(node1_reversed.encode("utf-8"))
    
    # Append the reversed sequence to node2
    node2.dna_bases = node2.dna_bases + node1_reversed_bits
    
    # Step 2: Update all paths
    for lineage, path in graph.paths.items():
        edges_to_remove = []
        
        for i, edge in enumerate(path.path_edges):
            # Remove the connecting edge between node1 and node2
            if edge.node1 == node1 and edge.node2 == node2:
                edges_to_remove.append(i)
                
            # Redirect ALL edges involving node1 to node2
            else:
                if edge.node1 == node1:
                    edge.node1 = node2
                    # Flip the side since we're now coming from the reversed part
                    edge.node1_side = not edge.node1_side
                if edge.node2 == node1:
                    edge.node2 = node2
                    # Flip the side since we're now going to the reversed part
                    edge.node2_side = not edge.node2_side
        
        # Remove connecting edges
        for i in reversed(edges_to_remove):
            del path.path_edges[i]
        
        path.node_count = len(path.path_edges)

def merge_nodes(graph):
    """
    Optimized merging pass on the graph with O(N + L×E) complexity per iteration.
    """
    initial_node_count = len(graph.nodes)
    print(f"Starting merge optimization on {initial_node_count} nodes")
    
    consumed_nodes = set()
    total_merges = 0
    iteration = 0
    
    while True:
        iteration += 1
        active_nodes = [n for n in graph.nodes if n not in consumed_nodes]
        print(f"\nIteration {iteration} : Active nodes: {len(active_nodes)}")
        
        # Step 1: Build adjacency information for this iteration
        print("Preprocessing graph")
        outgoing_from_true = defaultdict(set)   # node -> {(target_node, target_side)}
        incoming_to_side = defaultdict(set)     # (node, side) -> {(source_node, source_side)}
        
        edge_count = 0
        for lineage, path in graph.paths.items():
            for edge in path.path_edges:
                if edge.node1 not in consumed_nodes and edge.node2 not in consumed_nodes:
                    edge_count += 1
                    # Track outgoing edges from True side only
                    if edge.node1_side:  # From True side
                        outgoing_from_true[edge.node1].add((edge.node2, edge.node2_side))
                    
                    # Track all incoming edges
                    incoming_to_side[(edge.node2, edge.node2_side)].add((edge.node1, edge.node1_side))
        
        # Step 2: Find all valid merge candidates
        merge_candidates = []  # [(node1, node2, merge_type)]
        
        for node1 in active_nodes:
            # Get all targets from node1's True side
            targets = outgoing_from_true.get(node1, set())
            
            # Skip if no outgoing edges or multiple different targets
            if len(targets) != 1:
                continue
            
            target_node, target_side = list(targets)[0]
            
            # Skip if target already consumed
            if target_node in consumed_nodes:
                continue
            
            # Check if ALL edges to target's side come from node1's True side
            incoming = incoming_to_side.get((target_node, target_side), set())
            
            if incoming == {(node1, True)}:  # Only edges from node1's True side
                if target_side:
                    continue
                    merge_candidates.append((node1, target_node, "True-True"))
                else:
                    merge_candidates.append((node1, target_node, "True-False"))
        
        print(f"Processed {edge_count} edges, found {len(merge_candidates)} merge candidates")
        
        if not merge_candidates:
            print("No more merges possible")
            break
        
        # Step 3: Process merges (avoiding conflicts)
        print("Processing merges :")
        merges_this_iteration = 0
        nodes_modified_this_iteration = set()
        
        # Progress reporting setup
        total_candidates = len(merge_candidates)
        report_interval = max(1, total_candidates // 10)  # Report every 10%
        last_reported = 0
        
        for idx, (node1, node2, merge_type) in enumerate(merge_candidates):
            # Skip if either node was already involved in a merge this iteration
            if (node1 in consumed_nodes or node2 in consumed_nodes or 
                node1 in nodes_modified_this_iteration or node2 in nodes_modified_this_iteration):
                continue
            
            # Perform the merge
            if merge_type == "True-False":
                merge_node_True_False(node1, node2, graph)
            else:  # True-True
                merge_node_True_True(node1, node2, graph)
            
            consumed_nodes.add(node1)
            nodes_modified_this_iteration.add(node2)
            merges_this_iteration += 1
            total_merges += 1
            
            # Progress reporting
            if merges_this_iteration - last_reported >= report_interval:
                progress = ((idx + 1) / total_candidates) * 100
                print(f"  Progress: {merges_this_iteration} merges completed ({progress:.1f}%)")
                last_reported = merges_this_iteration
        
        # Final progress for this iteration if not already reported
        if merges_this_iteration > last_reported:
            print(f"  Progress: {merges_this_iteration} merges completed (100.0%)")
        
        # If we couldn't merge any candidates (all conflicted), we're done
        if merges_this_iteration == 0:
            print("No non-conflicting merges possible")
            break
    
    # Step 4: Final cleanup
    print("\nFinalizing...")
    graph.nodes -= consumed_nodes
    
    # Print summary
    final_node_count = len(graph.nodes)
    reduction = initial_node_count - final_node_count
    print(f"\n=== MERGE SUMMARY ===")
    print(f"Initial nodes: {initial_node_count}")
    print(f"Final nodes: {final_node_count}")
    print(f"Nodes merged: {reduction} ({reduction/initial_node_count*100:.1f}% reduction)")
    print(f"Total merges performed: {total_merges}")
    print(f"Total iterations: {iteration}")
    
    return (initial_node_count, final_node_count)