#!/usr/bin/env python3
"""
Simplified ARG Visualization Script

Creates specific visualizations for full ARGs without mutations.
Generates individual plots and consensus tree analysis.

Usage:
    python visualize_arg.py <input_trees_file> <output_folder>
"""

import sys
import os
import argparse
import tskit
import msprime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import pandas as pd
from pathlib import Path
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import defaultdict, Counter
from io_handler import MSerror, MSsuccess, MScompute

# Set matplotlib backend for non-interactive environments
plt.switch_backend('Agg')

# Consistent style settings
FIGURE_STYLE = {
    'figure.figsize': (12, 8),
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'font.weight': 'normal',
    'axes.titleweight': 'bold'
}
plt.rcParams.update(FIGURE_STYLE)


def get_node_types(ts):
    """Classify nodes into different types for full ARGs."""
    samples = set(ts.samples())
    re_nodes = set(nd.id for nd in ts.nodes() if nd.flags & msprime.NODE_IS_RE_EVENT)
    ca_nodes = set(nd.id for nd in ts.nodes() if nd.flags & msprime.NODE_IS_CA_EVENT)
    coalescent_nodes = set(range(ts.num_nodes)) - re_nodes - ca_nodes - samples
    
    return {
        'samples': samples,
        'recombination': re_nodes,
        'common_ancestor': ca_nodes,
        'coalescent': coalescent_nodes
    }


def create_local_trees_plot(ts, output_dir, basename):
    """Create local trees visualization with colored nodes."""
    MScompute("Creating local trees visualization...")
    
    # Get node types
    node_types = get_node_types(ts)
    
    # Get recombination times for grid lines
    re_times = [int(nd.time) for nd in ts.nodes() if nd.flags & msprime.NODE_IS_RE_EVENT]
    
    # Create style string for colored nodes with thicker edges
    style = ".edge {stroke-width: 3px} .y-axis .grid {stroke: #ff000033}"
    for u in node_types['recombination']:
        style += f".n{u} > .sym {{fill: red}}"
    for u in node_types['common_ancestor']:
        style += f".n{u} > .sym {{fill: orange}}"
    for u in node_types['coalescent']:
        style += f".n{u} > .sym {{fill: lightblue}}"
    
    # Create node labels for samples and recombination nodes
    node_labels = {u: str(u) for u in node_types['samples'] | node_types['recombination']}
    
    # Generate SVG
    svg_str = ts.draw_svg(
        size=(max(1200, ts.num_trees * 200), 500),
        y_axis=True,
        y_ticks=re_times if re_times else None,
        y_gridlines=bool(re_times),
        style=style,
        mutation_labels={},
        node_labels=node_labels
    )
    
    # Save SVG
    with open(os.path.join(output_dir, f"{basename}_local_trees.svg"), 'w') as f:
        f.write(svg_str)


def create_networkx_plot(ts, output_dir, basename):
    """Create NetworkX graph visualization."""
    MScompute("Creating NetworkX graph visualization...")
    
    # Convert to NetworkX graph
    D = dict(source=ts.edges_parent, target=ts.edges_child, 
             left=ts.edges_left, right=ts.edges_right)
    G = nx.from_pandas_edgelist(pd.DataFrame(D), edge_attr=True, 
                               create_using=nx.MultiDiGraph)
    
    nx.set_node_attributes(G, {n.id: {'flags': n.flags, 'time': n.time} 
                              for n in ts.nodes()})
    
    # Get node types for coloring
    node_types = get_node_types(ts)
    
    # Create layout using topological ordering
    for layer, nodes in enumerate(nx.topological_generations(G.reverse())):
        for node in nodes:
            G.nodes[node]["layer"] = layer
    
    pos = nx.multipartite_layout(G, subset_key="layer", align='horizontal')
    
    # Create figure
    plt.figure(figsize=(16, 10))
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.4, arrows=True, arrowsize=15, width=2)
    
    # Draw nodes with different colors
    if node_types['samples']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['samples']), 
                              node_color='#87CEEB', node_size=400, label='Samples')
    if node_types['recombination']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['recombination']), 
                              node_color='#FF4444', node_size=300, label='Recombination')
    if node_types['common_ancestor']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['common_ancestor']), 
                              node_color='#FFA500', node_size=250, label='Common Ancestor')
    if node_types['coalescent']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['coalescent']), 
                              node_color='#F08080', node_size=200, label='Coalescent')
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title(f"ARG Network Graph: {basename}", fontsize=16, fontweight='bold', pad=20)
    plt.legend(loc='upper right')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_networkx.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()


def create_tree_height_plot(ts, output_dir, basename):
    """Create tree height along genome plot."""
    MScompute("Creating tree height plot...")
    
    node_types = get_node_types(ts)
    
    plt.figure(figsize=(14, 6))
    
    positions = []
    heights = []
    for tree in ts.trees():
        positions.append(tree.interval.left)
        if tree.num_roots == 1:
            heights.append(ts.node(tree.root).time)
        else:
            heights.append(max(ts.node(root).time for root in tree.roots))
    
    plt.plot(positions, heights, 'b-', linewidth=3, label='Tree height', alpha=0.8)
    plt.fill_between(positions, heights, alpha=0.3, color='blue')
    
    # Add recombination events
    re_positions = []
    for re_node in node_types['recombination']:
        for tree in ts.trees():
            if re_node in tree.nodes():
                re_positions.append(tree.interval.left)
                break
    
    if re_positions:
        for pos in re_positions:
            plt.axvline(x=pos, color='red', linestyle='--', alpha=0.7, linewidth=2)
        plt.plot([], [], 'r--', label='Recombination events', linewidth=2)
    
    plt.title(f'Tree Height Along Genome: {basename}', fontsize=16, fontweight='bold')
    plt.xlabel('Genomic Position (bp)', fontsize=14)
    plt.ylabel('Height (generations)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_tree_height.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()

def compute_STEAC_matrix(ts):
    """
    Based on the Species Tree Estimation using Average Coalescence times method

    Compute average MRCA times between all pairs of samples in a tree sequence.

    This functions adjust the score based on the size of the locus.

    Parameters:
        ts (tskit.TreeSequence): A tree sequence with all genealogies.
    """
    samples = list(ts.samples())
    n_samples = len(samples)
    
    # Initialize an empty matrix to store coalescence times between each pair of samples
    coalescence_times = np.zeros((n_samples, n_samples))
    
    total_length = 0 

    for tree in ts.trees():
        span = tree.interval.right - tree.interval.left  # Locus lenght
        total_length += span  # Add to total genome length

        # Compare each pair of samples
        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):
                if i != j:  # Don't compare a sample to itself
                    # Find the MRCA of this pair in this tree
                    mrca = tree.mrca(sample1, sample2)
                    if mrca != tskit.NULL:
                        # Get the time of the MRCA
                        mrca_time = ts.node(mrca).time
                        # Add the time * span to the coalescence matrix
                        coalescence_times[i, j] += mrca_time * span

    # Average the coalescence times by dividing by the total genome length
    coalescence_times /= total_length

    return coalescence_times, samples


def create_hierarchical_clustering_plot(ts, output_dir, basename):
    """Create hierarchical clustering based on coalescence times."""
    MScompute("Creating hierarchical clustering plot...")
    
    coalescence_times, samples = compute_STEAC_matrix(ts)
    
    plt.figure(figsize=(12, 8))
    
    if len(samples) > 2:
        # Create condensed distance matrix (upper triangle)
        condensed_dist = []
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                condensed_dist.append(coalescence_times[i, j])
        
        if len(condensed_dist) > 0:
            # Perform hierarchical clustering
            linkage_matrix = linkage(condensed_dist, method='average')
            
            # Create dendrogram
            dendrogram(linkage_matrix, labels=[f'Sample {s}' for s in samples], 
                      leaf_font_size=12)
            plt.title(f'Hierarchical Clustering of Samples\n(Based on Coalescence Times): {basename}', 
                     fontsize=16, fontweight='bold')
            plt.ylabel('Coalescence Time', fontsize=14)
            plt.xlabel('Sample', fontsize=14)
    else:
        plt.text(0.5, 0.5, 'Need more than 2 samples for clustering', 
                transform=plt.gca().transAxes, ha='center', va='center', fontsize=14)
        plt.title(f'Hierarchical Clustering: {basename}', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_hierarchical_clustering.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="ARG visualization script")
    parser.add_argument("input_file", help="Input .trees file")
    parser.add_argument("output_dir", help="Output directory for visualizations")
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        MSerror(f"Input file {args.input_file} not found")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load tree sequence

    try:
        ts = tskit.load(args.input_file)
    except Exception as e:
        MSerror(f"Error loading tree sequence: {e}")
        sys.exit(1)
    
    # Get basename for output files
    basename = Path(args.input_file).stem
    
    
    # Create visualizations
    try:

        create_local_trees_plot(ts, args.output_dir, basename)
        create_networkx_plot(ts, args.output_dir, basename)
        create_tree_height_plot(ts, args.output_dir, basename)
        create_hierarchical_clustering_plot(ts, args.output_dir, basename)
        
        expected_files = [
            f"{basename}_local_trees.svg",
            f"{basename}_networkx.png",
            f"{basename}_tree_height.png",
            f"{basename}_hierarchical_clustering.png",
        ]
        
        for file in expected_files:
            if os.path.exists(os.path.join(args.output_dir, file)):
                pass
            else:
                MScompute(f"  ❌ {file} (not created)")
        
    except Exception as e:
        MSerror(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    MSsuccess("Arg visualizer")

if __name__ == "__main__":
    main()