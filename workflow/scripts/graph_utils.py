def merge_nodes(graph):
    """
    Merges nodes in the graph based on the following rule:
    - If a node has exactly one parent and that parent has exactly one child, merge them.
    
    Parameters:
        graph (Graph): The graph in which nodes will be merged.
    
    Returns:
        None (modifies the graph in place).
    """
    
    if not graph.nodes:
        raise ValueError("⚠️ MSpangepop -> Cannot merge nodes in an empty graph.")

    merged = set()
    
    # Traverse nodes in reverse order
    for node in reversed(graph.nodes):
        if node in merged:
            continue  # Skip already merged nodes
        
        if len(node.in_edges) == 1:
            parent = node.in_edges[0]
            if len(parent.out_edges) == 1:
                try:
                    parent.base += node.base
                    parent.out_edges = node.out_edges
                    for child in node.out_edges:
                        child.in_edges.remove(node)
                        child.in_edges.append(parent)

                    # Mark the merged node to avoid processing it again
                    merged.add(node)
                    graph.nodes.remove(node)

                except Exception as e:
                    raise RuntimeError(f"❌ MSpangepop -> Error merging nodes {parent.id} and {node.id}: {e}")

def save_to_gfa(graph, filename: str, sample: str, chromosome: str) -> None:
    """
    Saves the given graph to a GFA (Graphical Fragment Assembly) file.

    The GFA format consists of:
    - "S" lines defining sequence segments (nodes).
    - "L" lines defining links (edges) between segments.

    Parameters:
        graph (Graph): The graph to be saved.
        filename (str): The output GFA file path.

    Returns:
        None (writes to a file).
    """
    
    if not graph.nodes:
        raise ValueError("⚠️ MSpangepop -> Cannot save an empty graph to GFA.")

    try:
        with open(filename, 'w') as f:
            # Write nodes as segments in GFA format
            f.write(f"H\t{sample} chr_n{chromosome}\n")

            for node in graph.nodes:
                f.write(f"S\t{node.id}\t{node.base.decode()}\n")

            # Write edges (links) between nodes
            for node in graph.nodes:
                for out_node in node.out_edges:
                    f.write(f"L\t{node.id}\t+\t{out_node.id}\t+\t0M\n")

    except IOError as e:
        raise IOError(f"❌ MSpangepop -> Error writing to {filename}: {e}")
