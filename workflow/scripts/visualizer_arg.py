#!/usr/bin/env python3
"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Creates specific visualizations for full ARGs without mutations.
        Generates individual plots and consensus tree analysis.
"""

import sys
import os
import argparse
import tskit
import msprime
import numpy as np
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import networkx as nx
import pandas as pd
from pathlib import Path
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import defaultdict, Counter
from io_handler import MSerror, MSsuccess, MScompute, MSwarning
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from Bio.Phylo.Consensus import majority_consensus, strict_consensus
import io

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
    MScompute("Creating arg local trees visualization...")
    
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

def create_networkx_plot(ts, output_dir, basename, sample):
    """Create NetworkX graph visualization."""
    MScompute("Creating arg NetworkX graph visualization...")
    
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
    
    # Draw nodes with different colors - using darker colors
    if node_types['samples']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['samples']), 
                              node_color='#4682B4', node_size=400, label='Lineages')
    if node_types['recombination']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['recombination']), 
                              node_color='#8B0000', node_size=300, label='Recombination')
    if node_types['common_ancestor']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['common_ancestor']), 
                              node_color='#FF8C00', node_size=250, label='Common Ancestor')
    if node_types['coalescent']:
        nx.draw_networkx_nodes(G, pos, nodelist=list(node_types['coalescent']), 
                              node_color='#CD5C5C', node_size=200, label='Coalescent')
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=8)
    
    plt.title(f"ARG Network Graph: {sample} {basename}", fontsize=16, fontweight='bold', pad=20)
    plt.legend(loc='upper right')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_networkx.png"), 
                dpi=200, bbox_inches='tight')
    plt.close()


def create_tree_height_plot(ts, output_dir, basename, sample):
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
    
    plt.title(f'Tree Height Along Genome: {sample} {basename}', fontsize=16, fontweight='bold')
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
        span = tree.interval.right - tree.interval.left  # Locus length
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

def create_hierarchical_clustering_plot(ts, output_dir, basename, sample_name):
    """
    Create horizontal hierarchical clustering plot using STEAC with population-based coloring.
    Each branch is colored by the population of its descendants if all are the same.
    """
    MScompute("Creating STEAC horizontal hierarchical clustering plot...")

    coalescence_times, samples = compute_STEAC_matrix(ts)
    n = len(samples)

    if n > 2:
        # --- Get population information ---
        sample_populations = {}
        population_names = {}
        
        try:
            demography = ts.dump_tables().populations
            for i in range(len(demography)):
                if demography[i].metadata:
                    try:
                        import json
                        metadata = json.loads(demography[i].metadata.decode('utf-8'))
                        if 'name' in metadata:
                            population_names[i] = metadata['name']
                        else:
                            population_names[i] = f'Pop_{i}'
                    except:
                        population_names[i] = f'Pop_{i}'
                else:
                    population_names[i] = f'Pop_{i}'
        except:
            for i in range(ts.num_populations):
                population_names[i] = f'Pop_{i}'
        
        for sample_id in samples:
            pop_id = ts.node(sample_id).population
            sample_populations[sample_id] = pop_id
        
        # --- Color palette ---
        unique_pops = list(set(sample_populations.values()))
        dark_colors = ['#8B0000', '#00468B', '#2F4F4F', '#4B0082',
                       '#8B4513', '#191970', '#006400', '#8B008B']
        pop_colors = {pop: dark_colors[i % len(dark_colors)] for i, pop in enumerate(unique_pops)}
        
        # --- Labels ---
        labels = []
        label_colors = []
        sample_to_pop = []
        for s in samples:
            pop_id = sample_populations[s]
            pop_name = population_names.get(pop_id, f'Pop_{pop_id}')
            labels.append(f"L{s}")
            label_colors.append(pop_colors[pop_id])
            sample_to_pop.append(pop_id)
        
        # --- Distance matrix ---
        condensed_dist = [
            coalescence_times[i, j]
            for i in range(n)
            for j in range(i + 1, n)
        ]
        linkage_matrix = linkage(condensed_dist, method='average')
        
        # --- Build complete cluster population mapping ---
        def get_cluster_populations(cluster_id):
            """Get all populations in a cluster."""
            if cluster_id < n:
                # Leaf node
                return {sample_to_pop[cluster_id]}
            else:
                # Internal node
                idx = cluster_id - n
                left_child = int(linkage_matrix[idx, 0])
                right_child = int(linkage_matrix[idx, 1])
                left_pops = get_cluster_populations(left_child)
                right_pops = get_cluster_populations(right_child)
                return left_pops.union(right_pops)
        
        # Pre-compute all cluster populations
        cluster_pops = {}
        for i in range(n):
            cluster_pops[i] = {sample_to_pop[i]}
        for i in range(len(linkage_matrix)):
            cluster_pops[i + n] = get_cluster_populations(i + n)
        
        # --- Create figure ---
        fig, ax = plt.subplots(figsize=(max(12, n * 0.3), 10))
        
        # --- Plot dendrogram to get structure ---
        dend = dendrogram(
            linkage_matrix,
            labels=labels,
            leaf_rotation=0,
            leaf_font_size=8,
            orientation="left",
            ax=ax,
            color_threshold=0,  # All black by default
            above_threshold_color='black'
        )
        
        # Get the dendrogram structure
        icoord = dend['icoord']  # y-coordinates for horizontal dendrogram
        dcoord = dend['dcoord']  # x-coordinates (distances) for horizontal dendrogram
        
        # Get the correspondence between plot order and cluster IDs
        # The 'leaves' gives us the order of leaves in the plot
        # The color_list tells us which clusters are being plotted
        plot_order = dend['leaves']
        
        # Build a mapping from merge index to the corresponding icoord/dcoord index
        # This is needed because dendrogram may reorder things for visualization
        ivl = dend['ivl']  # The labels in plot order
        
        # Remove all existing lines
        for line in ax.get_lines():
            line.remove()
        
        # Create a mapping of which clusters appear at which y-coordinates
        # For leaves, we know their positions
        leaf_positions = {}
        y_tick_positions = [i * 10.0 for i in range(n)]  # Default spacing is 10
        for i, leaf_idx in enumerate(plot_order):
            leaf_positions[leaf_idx] = y_tick_positions[i]
        
        # Now redraw each segment with correct colors
        # We need to match each icoord/dcoord set with the correct merge
        for i in range(len(icoord)):
            y_coords = icoord[i]  # [left_child_y, merge_y, merge_y, right_child_y]
            x_coords = dcoord[i]  # [left_child_x, merge_x, merge_x, right_child_x]
            
            # Find which merge this corresponds to by looking at the structure
            # The merge that created this should be at distance x_coords[1] (or x_coords[2])
            merge_distance = x_coords[1]
            
            # Find the merge with this distance
            merge_idx = None
            for j in range(len(linkage_matrix)):
                if abs(linkage_matrix[j, 2] - merge_distance) < 1e-10:
                    # Check if the y-coordinates match what we expect
                    # This is the likely merge
                    merge_idx = j
                    break
            
            if merge_idx is not None:
                left_child = int(linkage_matrix[merge_idx, 0])
                right_child = int(linkage_matrix[merge_idx, 1])
                parent = merge_idx + n
                
                # Get populations for coloring
                left_pops = cluster_pops[left_child]
                right_pops = cluster_pops[right_child]
                parent_pops = cluster_pops[parent]
                
                # Determine colors - each segment is colored by what it leads TO
                left_color = pop_colors[next(iter(left_pops))] if len(left_pops) == 1 else 'black'
                right_color = pop_colors[next(iter(right_pops))] if len(right_pops) == 1 else 'black'
                parent_color = pop_colors[next(iter(parent_pops))] if len(parent_pops) == 1 else 'black'
            else:
                # Fallback to black if we can't match
                left_color = right_color = parent_color = 'black'
            
            # Draw the three segments
            # Left vertical branch (leads to left child)
            ax.plot([x_coords[0], x_coords[1]], [y_coords[0], y_coords[1]], 
                   color=left_color, lw=3)
            
            # Horizontal connector (represents the parent cluster)
            ax.plot([x_coords[1], x_coords[2]], [y_coords[1], y_coords[2]], 
                   color=parent_color, lw=3)
            
            # Right vertical branch (leads to right child)
            ax.plot([x_coords[2], x_coords[3]], [y_coords[2], y_coords[3]], 
                   color=right_color, lw=3)
        
        # Handle the root extension if there is one
        # Sometimes there's an additional vertical line extending from the root
        # This should be colored based on what the entire tree contains
        all_pops = set(sample_to_pop)
        root_color = pop_colors[next(iter(all_pops))] if len(all_pops) == 1 else 'black'
        
        # Check if there's a line at the far left that needs coloring
        max_x = max(max(d) for d in dcoord)
        for i in range(len(dcoord)):
            # Look for vertical segments at the maximum distance
            if abs(dcoord[i][0] - max_x) < 1e-10 or abs(dcoord[i][3] - max_x) < 1e-10:
                # This might be a root extension
                # Color it based on what it contains
                y1, y2 = min(icoord[i]), max(icoord[i])
                
                # Determine what this segment contains by its y-position
                # Find all leaves in this y-range
                contained_leaves = []
                for leaf_idx, y_pos in leaf_positions.items():
                    if y1 <= y_pos <= y2:
                        contained_leaves.append(leaf_idx)
                
                if contained_leaves:
                    contained_pops = {sample_to_pop[idx] for idx in contained_leaves}
                    segment_color = pop_colors[next(iter(contained_pops))] if len(contained_pops) == 1 else 'black'
                    
                    # Redraw this segment if needed
                    if dcoord[i][0] == max_x and dcoord[i][1] == max_x:
                        # Left side extension
                        ax.plot([dcoord[i][0], dcoord[i][1]], [icoord[i][0], icoord[i][1]], 
                               color=segment_color, lw=3)
        
        # Color the labels
        ylabels = ax.get_yticklabels()
        for i, label in enumerate(ylabels):
            label_text = label.get_text()
            for idx, orig_label in enumerate(labels):
                if label_text == orig_label:
                    label.set_color(label_colors[idx])
                    label.set_weight('bold')
                    break
        
        # --- Add legend ---
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=pop_colors[pop], label=population_names.get(pop, f'Pop_{pop}'))
            for pop in unique_pops
        ]
        ax.legend(handles=legend_elements, loc='upper left', title='Populations', frameon=True, 
                 fancybox=True, shadow=True)
        
        # --- Set titles and labels ---
        ax.set_title(
            f"STEAC Hierarchical Clustering\n(Species Tree Estimation using Average Coalescence times): {sample_name} {basename}",
            fontsize=14,
            fontweight='bold',
            pad=20
        )
        ax.set_xlabel("Average Coalescence Time", fontsize=12)
        ax.set_ylabel("Lineages", fontsize=12)
        
        # Add grid for better readability
        ax.grid(True, axis='x', alpha=0.3, linestyle='--')
        
        # Adjust layout
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{basename}_STEAC_hierarchical_clustering.png"),
                    dpi=1000, bbox_inches='tight')
        plt.close()
        
    else:
        plt.figure(figsize=(6, 4))
        plt.text(0.5, 0.5, 'Need more than 2 samples for clustering',
                 transform=plt.gca().transAxes, ha='center', va='center', fontsize=14)
        plt.title(f'STEAC Hierarchical Clustering: {sample_name} {basename}', fontsize=16, fontweight='bold')
        plt.savefig(os.path.join(output_dir, f"{basename}_STEAC_hierarchical_clustering.png"),
                    dpi=1000, bbox_inches='tight')
        plt.close()
        
def create_clustermap_with_dendrogram(ts, output_dir, basename, sample_name):
    """
    Create a clustered heatmap (clustermap) using STEAC coalescence times.
    Only shows left-side dendrogram; annotations are rounded to integers.
    Suitable for large sample sizes.
    """
    MScompute("Creating STEAC clustermap with left dendrogram...")

    coalescence_times, samples = compute_STEAC_matrix(ts)
    n = len(samples)

    if n > 2:
        # Use L instead of Sample for labels
        sample_labels = [f"L{s}" for s in samples]

        # Compute condensed distance matrix
        condensed_dist = [
            coalescence_times[i, j]
            for i in range(n)
            for j in range(i + 1, n)
        ]

        row_linkage = linkage(condensed_dist, method='average')

        df = pd.DataFrame(coalescence_times, index=sample_labels, columns=sample_labels)

        g = sns.clustermap(
            df,
            row_linkage=row_linkage,
            col_cluster=False,  # Only cluster rows
            cmap="coolwarm",
            figsize=(max(12, n * 0.35), max(10, n * 0.3)),
            annot=True,
            fmt=".0f",  # Round to integer
            annot_kws={"size": 6},
            dendrogram_ratio=(0.25, 0.01),  # Hide top dendrogram
            cbar_pos=(0.02, 0.8, 0.03, 0.18),
            linewidths=0.4,
            xticklabels=True,
            yticklabels=True,
            tree_kws={"linewidths": 1.5}
        )

        g.cax.set_visible(False)  # Hide colorbar

        plt.suptitle(f"STEAC Clustering (Left Tree Only): {sample_name} {basename}", fontsize=16, y=0.98)
        g.figure.subplots_adjust(top=0.93)
        g.savefig(os.path.join(output_dir, f"{basename}_STEAC_clustermap.png"),
                  dpi=1000, bbox_inches='tight')
        plt.close()
    else:
        MSwarning("Not enough samples for clustering heatmap")


def main():
    parser = argparse.ArgumentParser(description="ARG visualization script")
    parser.add_argument("input_file", help="Input .trees file")
    parser.add_argument("output_dir", help="Output directory for visualizations")
    parser.add_argument("chromosome")
    parser.add_argument("sample")
    args = parser.parse_args()
    sample = args.sample
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
    
    basename = f"chromosome_{args.chromosome}"
    
    try:
        create_local_trees_plot(ts, args.output_dir, basename)
        create_networkx_plot(ts, args.output_dir, basename, sample)
        create_tree_height_plot(ts, args.output_dir, basename, sample)
        
        # STEAC method - hierarchical clustering based on average coalescence times
        create_hierarchical_clustering_plot(ts, args.output_dir, basename, sample)
        create_clustermap_with_dendrogram(ts, args.output_dir, basename, sample)
        
        expected_files = [
            f"{basename}_local_trees.svg",
            f"{basename}_networkx.png",
            f"{basename}_tree_height.png",
            f"{basename}_STEAC_hierarchical_clustering.png",
            f"{basename}_STEAC_clustermap.png"
        ]
        
        for file in expected_files:
            full_path = os.path.join(args.output_dir, file)
            if os.path.exists(full_path):
                pass
            else:
                MSwarning(f"{file} (not created)")
    
    except Exception as e:
        MSerror(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    MSsuccess("Arg visualizations done!")


if __name__ == "__main__":
    main()