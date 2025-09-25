#!/usr/bin/env python3
"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage: Creates specific visualizations for full ARGs without mutations.
       Generates plots with optimized STEAC matrix computation.
"""

import sys
import os
import argparse
import time
import json
import tskit
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import msprime
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import stats
from collections import defaultdict
from io_handler import MSerror, MSsuccess, MScompute, MSwarning

# Set matplotlib backend for non-interactive environments
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
plt.switch_backend('Agg')

# Consistent color scheme - dark blue theme
COLOR_SCHEME = {
    'primary': '#1e3a5f',      # Dark blue
    'secondary': '#8B0000',     # Dark red
    'accent': '#FF8C00',        # Dark orange
    'light': '#4682B4',         # Steel blue
    'dark': '#0a1929',          # Very dark blue
    'grid': '#e0e0e0',          # Light gray for grids
    'background': '#f8f9fa'     # Light background
}

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

def compute_STEAC_matrix_original(ts):
    """
    O(T × n²) complexity.
    """
    samples = list(ts.samples())
    n_samples = len(samples)
    coalescence_times = np.zeros((n_samples, n_samples))
    total_length = 0

    for tree in ts.trees():
        span = tree.interval.right - tree.interval.left
        total_length += span

        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):
                if i != j:
                    mrca = tree.mrca(sample1, sample2)
                    if mrca != tskit.NULL:
                        mrca_time = ts.node(mrca).time
                        coalescence_times[i, j] += mrca_time * span

    coalescence_times /= total_length
    return coalescence_times, samples

def compute_STEAC_matrix_vectorized(ts):
    """
    Vectorized STEAC computation using numpy operations.
    More efficient for larger sample sizes.
    """
    samples = np.array(ts.samples())
    n_samples = len(samples)
    coalescence_times = np.zeros((n_samples, n_samples))
    
    for tree in ts.trees():
        span = tree.interval.right - tree.interval.left
        
        # Pre-compute all MRCAs for this tree
        mrca_matrix = np.zeros((n_samples, n_samples), dtype=np.int32)
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                mrca = tree.mrca(samples[i], samples[j])
                if mrca != tskit.NULL:
                    mrca_matrix[i, j] = mrca_matrix[j, i] = mrca
        
        # Vectorized time lookup and accumulation
        mask = mrca_matrix != 0
        times = np.zeros_like(mrca_matrix, dtype=np.float64)
        times[mask] = ts.nodes_time[mrca_matrix[mask]]
        coalescence_times += times * span
    
    coalescence_times /= ts.sequence_length
    return coalescence_times, list(samples)

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


def create_tree_height_plot(ts, output_dir, basename, sample):
    """Create tree height along genome plot as histogram bars."""
    MScompute("Creating tree height plot...")
    
    plt.figure(figsize=(14, 6))
    
    # Collect data for bars AND boundaries in single pass
    left_positions = []
    widths = []
    heights = []
    recomb_positions = set()
    
    for i, tree in enumerate(ts.trees()):
        left_positions.append(tree.interval.left)
        widths.append(tree.interval.right - tree.interval.left)
        
        if tree.num_roots == 1:
            heights.append(ts.node(tree.root).time)
        else:
            heights.append(max(ts.node(root).time for root in tree.roots))
        
        # Add tree boundaries (except for first and last)
        if i > 0:  # Not the first tree
            recomb_positions.add(tree.interval.left)
    
    # Create bar plot with dark blue color
    plt.bar(left_positions, heights, width=widths, 
            align='edge',
            color='#1e3a5f',
            edgecolor='#0a1929',
            linewidth=0.3,
            alpha=0.9,
            label='Tree height')
    
    # Add recombination breakpoints
    if recomb_positions:
        for pos in sorted(recomb_positions):
            plt.axvline(x=pos, color='#8b0000', linestyle='--', alpha=0.8, linewidth=.5)
        plt.plot([], [], color='#8b0000', linestyle='--', linewidth=2, label='Recombination breakpoints')
    
    plt.title(f'Tree Height Along Genome: {sample} {basename}', fontsize=16, fontweight='bold')
    plt.xlabel('Genomic Position (bp)', fontsize=14)
    plt.ylabel('Height (generations)', fontsize=14)
    plt.legend(loc='best', frameon=True, fancybox=True, shadow=True)
    plt.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    plt.ylim(bottom=0)
    
    ax = plt.gca()
    ax.set_facecolor('#f8f9fa')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_tree_height.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def create_hierarchical_clustering_plot(ts, output_dir, basename, sample_name, coalescence_times, samples):
    """
    Create horizontal hierarchical clustering plot using STEAC with population-based coloring.
    Each branch is colored by the population of its descendants if all are the same.
    """
    MScompute("Creating STEAC horizontal hierarchical clustering plot...")

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

def create_clustermap_with_dendrogram(ts, output_dir, basename, sample_name, coalescence_times, samples):
    """
    Create a clustered heatmap using STEAC coalescence times.
    """
    MScompute("Creating STEAC clustermap...")

    n = len(samples)

    if n > 2:
        sample_labels = [f"L{s}" for s in samples]

        # Compute condensed distance matrix
        condensed_dist = [
            coalescence_times[i, j]
            for i in range(n)
            for j in range(i + 1, n)
        ]

        row_linkage = linkage(condensed_dist, method='average')
        df = pd.DataFrame(coalescence_times, index=sample_labels, columns=sample_labels)

        # Create custom colormap using our color scheme
        cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
            'custom', 
            ['white', COLOR_SCHEME['light'], COLOR_SCHEME['primary'], COLOR_SCHEME['dark']]
        )

        g = sns.clustermap(
            df,
            row_linkage=row_linkage,
            col_cluster=False,
            cmap=cmap,
            figsize=(max(12, n * 0.35), max(10, n * 0.3)),
            annot=True,
            fmt=".0f",
            annot_kws={"size": 6},
            dendrogram_ratio=(0.25, 0.01),
            cbar_pos=(0.02, 0.8, 0.03, 0.18),
            linewidths=0.4,
            xticklabels=True,
            yticklabels=True,
            tree_kws={"linewidths": 1.5, "colors": COLOR_SCHEME['dark']}
        )

        g.cax.set_visible(False)
        
        plt.suptitle(f"STEAC Clustering: {sample_name} {basename}", 
                     fontsize=16, fontweight='bold', y=0.98)
        g.figure.subplots_adjust(top=0.93)
        g.savefig(os.path.join(output_dir, f"{basename}_STEAC_clustermap.png"),
                  dpi=300, bbox_inches='tight')
        plt.close()
    else:
        MSwarning("Not enough samples for clustering heatmap")


def plot_tree_age_vs_width(ts, output_dir, basename, sample):
    """Plot tree age versus tree width."""
    MScompute("Creating tree age vs width plot...")
    
    # Collect data
    widths = []
    ages = []
    
    for tree in ts.trees():
        width = tree.interval.right - tree.interval.left
        widths.append(width)
        
        if tree.num_roots == 1:
            age = ts.node(tree.root).time
        else:
            age = max(ts.node(root).time for root in tree.roots)
        ages.append(age)
    
    widths = np.array(widths)
    ages = np.array(ages)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Scatter plot
    scatter = ax1.scatter(widths, ages,
                         alpha=0.6,
                         s=50,
                         c=ages,
                         cmap='viridis',
                         edgecolors=COLOR_SCHEME['dark'],
                         linewidth=0.5)
    
    cbar = plt.colorbar(scatter, ax=ax1)
    cbar.set_label('Tree Age (generations)', fontsize=11)
    
    # Calculate correlation
    correlation, p_value = stats.pearsonr(widths, ages)
    ax1.text(0.05, 0.95, f'Pearson r = {correlation:.3f}\np-value = {p_value:.2e}',
             transform=ax1.transAxes,
             bbox=dict(boxstyle='round', facecolor=COLOR_SCHEME['background'], 
                      edgecolor=COLOR_SCHEME['primary'], alpha=0.9),
             verticalalignment='top',
             fontsize=10)
    
    # Add trend line
    z = np.polyfit(widths, ages, 1)
    p = np.poly1d(z)
    x_trend = np.linspace(widths.min(), widths.max(), 100)
    ax1.plot(x_trend, p(x_trend), color=COLOR_SCHEME['secondary'], 
             linestyle='--', alpha=0.7, linewidth=2, label='Linear trend')
    
    ax1.set_xlabel('Tree Width (bp)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Tree Age (generations)', fontsize=12, fontweight='bold')
    ax1.set_title('Tree Age vs Genomic Width', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':', color=COLOR_SCHEME['grid'])
    ax1.set_facecolor(COLOR_SCHEME['background'])
    ax1.legend()
    
    # Hexbin density plot with custom colormap
    cmap_hex = plt.cm.colors.LinearSegmentedColormap.from_list(
        'custom_hex', 
        ['white', COLOR_SCHEME['light'], COLOR_SCHEME['accent'], COLOR_SCHEME['secondary']]
    )
    hexbin = ax2.hexbin(widths, ages, gridsize=30, cmap=cmap_hex, mincnt=1)
    
    cbar2 = plt.colorbar(hexbin, ax=ax2)
    cbar2.set_label('Count', fontsize=11)
    
    ax2.set_xlabel('Tree Width (bp)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Tree Age (generations)', fontsize=12, fontweight='bold')
    ax2.set_title('Density Plot: Tree Age vs Width', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':', color=COLOR_SCHEME['grid'])
    ax2.set_facecolor(COLOR_SCHEME['background'])
    
    fig.suptitle(f'Tree Age-Width Relationship: {sample} {basename}',
                 fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_tree_age_vs_width.png"),
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # Print summary statistics
    MScompute(f"Tree statistics:")
    MScompute(f"  - Number of trees: {len(widths)}")
    MScompute(f"  - Mean tree width: {np.mean(widths):.2f} bp")
    MScompute(f"  - Mean tree age: {np.mean(ages):.2f} generations")
    MScompute(f"  - Correlation (age vs width): {correlation:.3f}")


def main():
    parser = argparse.ArgumentParser(description="ARG visualization script with vectorized STEAC computation")
    parser.add_argument("input_file", help="Input .trees file")
    parser.add_argument("output_dir", help="Output directory for visualizations")
    parser.add_argument("chromosome", help="Chromosome identifier")
    parser.add_argument("sample", help="Sample identifier")
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
        MScompute(f"Loaded tree sequence: {ts.num_trees} trees, {ts.num_samples} samples")
    except Exception as e:
        MSerror(f"Error loading tree sequence: {e}")
        sys.exit(1)

    basename = f"chromosome_{args.chromosome}"
    sample = args.sample

    try:
        # Compute STEAC matrix using vectorized method
        MScompute("Computing STEAC matrix")
        coalescence_times, samples = compute_STEAC_matrix_vectorized(ts)

        # Create all original visualizations
        create_tree_height_plot(ts, args.output_dir, basename, sample)
        plot_tree_age_vs_width(ts, args.output_dir, basename, sample)
        create_hierarchical_clustering_plot(ts, args.output_dir, basename, sample, 
                                            coalescence_times, samples)
        create_clustermap_with_dendrogram(ts, args.output_dir, basename, sample, 
                                          coalescence_times, samples)

    except Exception as e:
        MSerror(f"Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    MSsuccess("All ARG visualizations complete!")


if __name__ == "__main__":
    main()