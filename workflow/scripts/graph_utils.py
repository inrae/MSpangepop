import random
from io_handler import MSerror


def mutate_base(original_base: str, traition_matrix) -> str:
    """Uses the provided transition matrix to determine the mutated base."""
    return random.choices(
        population=list(traition_matrix[original_base].keys()),
        weights=list(traition_matrix[original_base].values()),
        k=1
    )[0]

def generate_sequence(length, transition_matrix, start_frequencies=None):
    """
    Generate a DNA sequence of given length using a transition matrix.
    
    Args:
        length (int): Desired sequence length
        transition_matrix (dict): Nested dict with transition probabilities
                                 e.g., {'A': {'A': 0.3, 'T': 0.2, 'G': 0.3, 'C': 0.2}, ...}
        start_frequencies (dict, optional): Initial nucleotide probabilities
                                          e.g., {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}
                                          If None, uses uniform distribution
    
    Returns:
        str: Generated DNA sequence
    """
    if length <= 0:
        return ""
    
    # Default to uniform start frequencies if not provided
    if start_frequencies is None:
        nucleotides = list(transition_matrix.keys())
        start_frequencies = {base: 1.0/len(nucleotides) for base in nucleotides}
    
    # Choose starting nucleotide
    current_base = random.choices(
        list(start_frequencies.keys()),
        weights=list(start_frequencies.values())
    )[0]
    
    sequence = [current_base]
    
    # Generate subsequent bases using transition matrix
    for _ in range(length - 1):
        next_base = random.choices(
            list(transition_matrix[current_base].keys()),
            weights=list(transition_matrix[current_base].values())
        )[0]
        sequence.append(next_base)
        current_base = next_base
    
    return ''.join(sequence)

'''
def merge_nodes(graph):
    """
    Merges nodes in the graph based on the following rule:
    - If a node has exactly one parent and that parent has exactly one child, merge them.
    - Updates paths accordingly to reflect merged nodes.

    Parameters:
        graph (Graph): The graph in which nodes will be merged.

    Returns:
        None (modifies the graph in place).
    """
    
    if not graph.nodes:
        raise MSerror("Cannot merge nodes in an empty graph.")

    merged_nodes = set()  # Track nodes to remove

    for node in reversed(graph.nodes):
        if node in merged_nodes:
            continue  # Skip already merged nodes

        if len(node.in_edges) == 1:
            parent = node.in_edges[0]
            if len(parent.out_edges) == 1:
                # Merge node into its parent
                parent.base += node.base
                parent.out_edges = node.out_edges
                
                for child in node.out_edges:
                    child.in_edges.remove(node)
                    child.in_edges.append(parent)

                # Update paths: replace merged node with its parent
                for path in graph.paths.values():
                    if node in path.nodes:
                        node_index = path.nodes.index(node)
                        path.nodes[node_index] = parent  # Replace with merged parent
                        path.nodes = list(dict.fromkeys(path.nodes))  # Remove duplicates

                # Mark node for removal
                merged_nodes.add(node)

    # Remove merged nodes in one go (avoids O(n²) slow removals)
    graph.nodes = [node for node in graph.nodes if node not in merged_nodes]

def save_to_gfa(graph, filename: str, sample: str, chromosome: str) -> None:
    """
    Saves the given graph to a GFA (Graphical Fragment Assembly) file.

    The GFA format consists of:
    - "S" lines defining sequence segments (nodes).
    - "L" lines defining links (edges) between segments.
    - "P" lines defining paths.

    Parameters:
        graph (Graph): The graph to be saved.
        filename (str): The output GFA file path.

    Returns:
        None (writes to a file).
    """
    
    if not graph.nodes:
        raise MSerror("Cannot save an empty graph to GFA.")

    try:
        with open(filename, 'w') as f:
            # Write header
            f.write(f"H\t{sample} chr_n{chromosome}\n")

            # Write nodes (segments)
            for node in graph.nodes:
                f.write(f"S\t{node.id}\t{node}\n")

            # Write edges (links)
            for node in graph.nodes:
                for out_node in node.out_edges:
                    f.write(f"L\t{node.id}\t+\t{out_node.id}\t+\t0M\n")

            # Write paths in the format "P\t{path.lineage}\t{node1>node2>...}"
            for path in graph.paths.values():
                node_order = ">".join(str(node.id) for node in path.nodes)
                f.write(f"P\tlineage_{path.lineage}\t{node_order}\n")

    except IOError as e:
        raise MSerror(f"Error writing to {filename}: {e}")
'''