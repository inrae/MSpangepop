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
from io_handler import MSerror, MSsuccess, MScompute, MSwarning

# BioPython imports for consensus trees
try:
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Tree, Clade
    from Bio.Phylo.Consensus import majority_consensus, strict_consensus
    import io
except ImportError:
    MSerror("BioPython not installed. Please install with: pip install biopython")
    sys.exit(1)

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


def tskit_to_biopython_tree(tree, ts):
    """
    Convert a tskit tree to a BioPython tree.
    """
    def build_clade(node):
        """Recursively build Bio.Phylo clades from tskit nodes."""
        clade = Clade()
        
        # If it's a sample (leaf), give it a name
        if tree.is_sample(node):
            clade.name = f"Sample_{node}"
        
        # Get children
        children = tree.children(node)
        if children:
            for child in children:
                child_clade = build_clade(child)
                # Set branch length as the time difference
                child_clade.branch_length = ts.node(node).time - ts.node(child).time
                clade.clades.append(child_clade)
        
        return clade
    
    # Start from the root
    if tree.num_roots == 1:
        root_clade = build_clade(tree.root)
        bio_tree = Tree(root=root_clade)
        return bio_tree
    else:
        # Handle multiple roots by creating a super-root
        super_root = Clade()
        for root in tree.roots:
            root_clade = build_clade(root)
            root_clade.branch_length = 0
            super_root.clades.append(root_clade)
        bio_tree = Tree(root=super_root)
        return bio_tree


def create_majority_rule_consensus_plot(ts, output_dir, basename):
    """
    Create a majority-rule consensus tree visualization using BioPython.
    
    This method converts all local trees to BioPython format and uses
    the majority_consensus function to find the consensus tree.
    """
    MScompute("Creating majority-rule consensus tree using BioPython...")
    
    samples = list(ts.samples())
    n_samples = len(samples)
    
    # Convert all local trees to BioPython format
    bio_trees = []
    tree_weights = []
    
    for tree in ts.trees():
        # Weight by the span of the tree
        span = tree.interval.right - tree.interval.left
        weight = span / ts.sequence_length
        
        # Convert to BioPython tree
        bio_tree = tskit_to_biopython_tree(tree, ts)
        bio_trees.append(bio_tree)
        tree_weights.append(weight)
    
    # Create consensus tree (majority rule with 0.5 cutoff)
    try:
        # BioPython's majority_consensus expects equal weight trees,
        # so we'll duplicate trees based on their relative weights
        # Convert weights to approximate integer counts
        min_weight = min(tree_weights)
        scaled_counts = [int(round(w / min_weight)) for w in tree_weights]
        
        # Create weighted tree list
        weighted_trees = []
        for tree, count in zip(bio_trees, scaled_counts):
            weighted_trees.extend([tree] * count)
        
        # Calculate majority consensus (default cutoff is 0.5)
        consensus_tree = majority_consensus(weighted_trees, cutoff=0.5)
        
        # Also try strict consensus for comparison
        strict_tree = strict_consensus(bio_trees)
        
    except Exception as e:
        MScompute(f"  ⚠️ Error creating consensus tree: {e}")
        consensus_tree = None
        strict_tree = None
    
    # Create visualization
    fig = plt.figure(figsize=(16, 8))
    
    if consensus_tree is not None:
        # Plot majority-rule consensus
        ax1 = fig.add_subplot(1, 2, 1)
        try:
            Phylo.draw(consensus_tree, axes=ax1, do_show=False)
            ax1.set_title(f'Majority-Rule Consensus Tree (>50%)\n{basename}', 
                         fontsize=14, fontweight='bold')
        except Exception as e:
            ax1.text(0.5, 0.5, f'Could not draw majority consensus tree:\n{e}', 
                    transform=ax1.transAxes, ha='center', va='center')
            ax1.set_title('Majority-Rule Consensus Tree', fontsize=14, fontweight='bold')
        
        # Plot strict consensus for comparison
        ax2 = fig.add_subplot(1, 2, 2)
        if strict_tree is not None:
            try:
                Phylo.draw(strict_tree, axes=ax2, do_show=False)
                ax2.set_title(f'Strict Consensus Tree (100%)\n{basename}', 
                             fontsize=14, fontweight='bold')
            except:
                ax2.text(0.5, 0.5, 'Strict consensus tree is empty\n(no bipartitions in all trees)', 
                        transform=ax2.transAxes, ha='center', va='center')
                ax2.set_title('Strict Consensus Tree', fontsize=14, fontweight='bold')
        else:
            ax2.text(0.5, 0.5, 'No strict consensus tree\n(no common structure)', 
                    transform=ax2.transAxes, ha='center', va='center')
            ax2.set_title('Strict Consensus Tree', fontsize=14, fontweight='bold')
    else:
        plt.text(0.5, 0.5, 'Could not create consensus trees from ARG', 
                transform=plt.gca().transAxes, ha='center', va='center', fontsize=14)
    
    plt.suptitle(f'Consensus Trees from ARG: {basename}', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_majority_rule_consensus.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a separate detailed plot of just the majority consensus
    if consensus_tree is not None:
        plt.figure(figsize=(12, 10))
        try:
            # Draw with more detail
            Phylo.draw(consensus_tree, do_show=False, 
                      branch_labels=lambda c: f"{c.confidence:.0%}" if c.confidence else "")
            plt.title(f'Majority-Rule Consensus Tree (>50% bipartitions): {basename}', 
                     fontsize=16, fontweight='bold')
            plt.xlabel('Branch Length', fontsize=12)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{basename}_majority_consensus_detailed.png"), 
                        dpi=300, bbox_inches='tight')
        except:
            pass
        plt.close()
    
    # Also create a text summary
    with open(os.path.join(output_dir, f"{basename}_majority_rule_summary.txt"), 'w') as f:
        f.write(f"Consensus Tree Summary for {basename}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total number of local trees: {ts.num_trees}\n")
        f.write(f"Total sequence length: {ts.sequence_length}\n")
        f.write(f"Number of samples: {n_samples}\n\n")
        
        if consensus_tree is not None:
            f.write("Majority-Rule Consensus Tree (>50% support):\n")
            f.write("-" * 40 + "\n")
            # Write tree in Newick format
            newick_str = io.StringIO()
            Phylo.write(consensus_tree, newick_str, 'newick')
            f.write(f"Newick format: {newick_str.getvalue()}\n\n")
            
            # Count clades
            all_clades = list(consensus_tree.find_clades())
            terminal_clades = [c for c in all_clades if c.is_terminal()]
            internal_clades = [c for c in all_clades if not c.is_terminal() and c != consensus_tree.root]
            
            f.write(f"Number of clades: {len(all_clades)}\n")
            f.write(f"Number of terminal clades (samples): {len(terminal_clades)}\n")
            f.write(f"Number of internal clades: {len(internal_clades)}\n")
        else:
            f.write("Could not create consensus tree from ARG.\n")
            f.write("This may indicate very high recombination or incompatible tree structures.\n")


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


def create_hierarchical_clustering_plot(ts, output_dir, basename):
    """
    Create hierarchical clustering based on STEAC (Species Tree Estimation using Average Coalescence times).
    
    This method uses the average coalescence times between samples across all local trees
    to build a distance matrix and perform hierarchical clustering.
    """
    MScompute("Creating STEAC hierarchical clustering plot...")
    
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
            plt.title(f'STEAC Hierarchical Clustering\n(Species Tree Estimation using Average Coalescence times): {basename}', 
                     fontsize=16, fontweight='bold')
            plt.ylabel('Average Coalescence Time', fontsize=14)
            plt.xlabel('Sample', fontsize=14)
    else:
        plt.text(0.5, 0.5, 'Need more than 2 samples for clustering', 
                transform=plt.gca().transAxes, ha='center', va='center', fontsize=14)
        plt.title(f'STEAC Hierarchical Clustering: {basename}', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{basename}_STEAC_hierarchical_clustering.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()


def create_clustermap_with_dendrogram(ts, output_dir, basename):
    """Create a clustered heatmap with dendrogram using STEAC method."""
    MScompute("Creating STEAC clustered heatmap with dendrogram...")

    coalescence_times, samples = compute_STEAC_matrix(ts)
    sample_labels = [f"Sample {s}" for s in samples]

    if len(samples) > 2:
        # Convert to condensed distance matrix for clustering
        condensed_dist = []
        for i in range(len(samples)):
            for j in range(i+1, len(samples)):
                condensed_dist.append(coalescence_times[i, j])
        
        # Pre-compute the linkage
        row_linkage = linkage(condensed_dist, method='average')
        col_linkage = row_linkage  # Use same linkage for rows and columns (symmetric matrix)
        
        # Create DataFrame for heatmap
        df = pd.DataFrame(coalescence_times, index=sample_labels, columns=sample_labels)

        # Create clustermap with pre-computed linkage
        g = sns.clustermap(df,
                           row_linkage=row_linkage,
                           col_linkage=col_linkage,
                           cmap='coolwarm',
                           figsize=(14, 12),
                           annot=True,
                           fmt=".2f",
                           dendrogram_ratio=(0.25, 0.25),
                           cbar_pos=(0.02, 0.8, 0.03, 0.18),
                           linewidths=0.5,
                           tree_kws={'linewidths': 3.0},  # Increased tree branch width
                           )
        
        # Remove the legend/colorbar
        g.cax.set_visible(False)
        
        plt.suptitle(f"STEAC Clustering Analysis: {basename}", fontsize=18, y=0.98)
        g.fig.subplots_adjust(top=0.94)
        g.savefig(os.path.join(output_dir, f"{basename}_STEAC_clustermap.png"),
                  dpi=300, bbox_inches='tight')
        plt.close()
    else:
        MScompute("  ⚠️ Not enough samples for clustering heatmap")


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
    
    basename = Path(args.input_file).stem
    
    try:
        create_local_trees_plot(ts, args.output_dir, basename)
        create_networkx_plot(ts, args.output_dir, basename)
        create_tree_height_plot(ts, args.output_dir, basename)
        
        # STEAC method - hierarchical clustering based on average coalescence times
        create_hierarchical_clustering_plot(ts, args.output_dir, basename)
        create_clustermap_with_dendrogram(ts, args.output_dir, basename)
        
        # Majority-rule consensus tree method
        #create_majority_rule_consensus_plot(ts, args.output_dir, basename)

        expected_files = [
            f"{basename}_local_trees.svg",
            f"{basename}_networkx.png",
            f"{basename}_tree_height.png",
            f"{basename}_STEAC_hierarchical_clustering.png",
            f"{basename}_STEAC_clustermap.png"
            #f"{basename}_majority_rule_consensus.png",
            #f"{basename}_majority_consensus_detailed.png",
            #f"{basename}_majority_rule_summary.txt",
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