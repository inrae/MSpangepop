"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

This scripts holds usefull function for the graph creation step 
"""

import random
from io_handler import MSerror, MSsuccess, MSwarning, MScompute
from datetime import datetime
import os
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
import matplotlib.pyplot as plt # type: ignore
import matplotlib.patches as mpatches # type: ignore
import numpy as np
from collections import defaultdict
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

def mutate_base(original_base: str, traition_matrix: dict) -> str:
    """Uses the provided transition matrix to determine the mutated base."""
    return random.choices(
        population=list(traition_matrix[original_base].keys()),
        weights=list(traition_matrix[original_base].values()),
        k=1
    )[0]

def generate_sequence(length: int, transition_matrix: dict) -> list:
    """
    Generate a DNA sequence of given length using a transition matrix.
    """
    if length <= 0:
        return []
    
    start_frequencies = transition_matrix["_"] # The start frequencie is provided as "_"
    
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
    
    return sequence

def gather_lineages(traversal):
    """Extract lineages for each tree and return as a list to maintain order."""
    tree_lineages = []
    for tree in traversal:
        tree_index = tree.get("tree_index", "unknown")
        lineages = set(tree.get("lineages", []))  # Convert to set for build_from_sequence
        tree_lineages.append((tree_index, lineages))
    return tree_lineages

class MutationRecap:
    """Tracks all mutation applications for recap file generation."""
    
    def __init__(self, sample: str, chromosome: str):
        self.sample = sample
        self.chromosome = chromosome
        self.start_time = datetime.now()
        self.mutations = []
        self.summary = {
            "total_attempted": 0,
            "total_successful": 0,
            "total_failed": 0,
            "by_type": {}
        }
    
    def add_mutation(self, tree_index: int, node_id: int, mutation_type: str, 
                     position: int, length: int, affected_lineages: set, 
                     success: bool, error_msg: str = None):
        """Record a mutation attempt."""
        self.mutations.append({
            "tree_index": tree_index,
            "node_id": node_id,
            "type": mutation_type,
            "position": position,
            "length": length,
            "affected_lineages": sorted(list(affected_lineages)),
            "success": success,
            "error": error_msg
        })
        
        # Update summary
        self.summary["total_attempted"] += 1
        if success:
            self.summary["total_successful"] += 1
        else:
            self.summary["total_failed"] += 1
        
        # Track by type
        if mutation_type not in self.summary["by_type"]:
            self.summary["by_type"][mutation_type] = {"attempted": 0, "successful": 0, "failed": 0}
        
        self.summary["by_type"][mutation_type]["attempted"] += 1
        if success:
            self.summary["by_type"][mutation_type]["successful"] += 1
        else:
            self.summary["by_type"][mutation_type]["failed"] += 1
    
    def save_recap(self, filepath: str):
        """Save recap to file with inline format."""
        with open(filepath, 'w') as f:
            # Write header
            f.write("🔹 MSpangepop Mutation Recap File\n")
            f.write("This file recaps all mutations added in the graph.\n")
            f.write("-"*120+"\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Sample: {self.sample}\n")
            f.write(f"Chromosome: {self.chromosome}\n")
            f.write(f"Processing Duration: {datetime.now() - self.start_time}\n")
            f.write("-"*120+"\n")
            
            # Write summary
            f.write("SUMMARY\n")
            f.write(f"Total Mutations Attempted: {self.summary['total_attempted']}\n")
            f.write(f"Total Successful: {self.summary['total_successful']}\n")
            f.write(f"Total Failed: {self.summary['total_failed']}\n")
            if self.summary['total_attempted'] > 0:
                success_rate = self.summary['total_successful'] / self.summary['total_attempted'] * 100
                f.write(f"Success Rate: {success_rate:.2f}%\n\n")
            else:
                f.write("Success Rate: N/A\n\n")
            
            f.write("-"*120+"\n")
            # Summary by type
            f.write("By Mutation Type:\n")
            for mut_type, stats in sorted(self.summary["by_type"].items()):
                success_rate = stats['successful'] / max(1, stats['attempted']) * 100
                f.write(f"{mut_type:8} - Attempted: {stats['attempted']:4d}, "
                       f"Successful: {stats['successful']:4d}, "
                       f"Failed: {stats['failed']:4d} "
                       f"(Success Rate: {success_rate:6.2f}%)\n")
            f.write("\n")
            f.write("-"*120+"\n")
            # Write column headers for mutations
            f.write("DETAILED MUTATION LOG\n")
            f.write("-" * 120 + "\n")
            f.write(f"{'Mutation #':<12}\t{'Tree Index':<12}\t{'Tree Node':<12}\t"
                   f"{'Type':<8}\t{'Position':<10}\t{'Length':<8}\t"
                   f"{'Lineages':<30}\t{'Status':<10}\n")
            f.write("-" * 120 + "\n")
            
            # Write mutations in inline format
            for i, mut in enumerate(self.mutations, 1):
                lineages_str = ','.join(map(str, mut['affected_lineages']))
                if len(lineages_str) > 30:
                    lineages_str = lineages_str[:27] + "..."
                
                # Handle None values properly
                tree_index_str = str(mut['tree_index']) if mut['tree_index'] is not None else "N/A"
                node_id_str = str(mut['node_id']) if mut['node_id'] is not None else "N/A"
                type_str = str(mut['type']) if mut['type'] is not None else "N/A"
                position_str = str(mut['position']) if mut['position'] is not None else "N/A"
                length_str = str(mut['length']) if mut['length'] is not None else "N/A"
                status = 'SUCCESS' if mut['success'] else 'FAILED'
                
                f.write(f"{i:<12}\t{tree_index_str:<12}\t{node_id_str:<12}\t"
                       f"{type_str:<8}\t{position_str:<10}\t{length_str:<8}\t"
                       f"{lineages_str:<30}\t{status:<10}")
                
                # Add error message on the same line if failed
                if mut['error'] and not mut['success']:
                    f.write(f"\t# {mut['error']}")
                
                f.write("\n")
            
            f.write("-" * 120 + "\n")
            f.write(f"End of recap - Total mutations: {len(self.mutations)}\n")

class VariantSizeVisualizer:
    """Creates visualization graphs for variant sizes from mutation data."""
    
    def __init__(self, sample: str, chromosome: str, reference_length: int = None):
        self.sample = sample
        self.chromosome = chromosome
        self.reference_length = reference_length  # Original/reference genome length
        self.variant_sizes = defaultdict(list)
        self.variant_positions = defaultdict(list)
        self.lineage_lengths = {}
        self.max_position = 0
        # Track variants by lineage for shared variant analysis
        self.variants_by_lineage = defaultdict(set)
        self.variant_lineages = defaultdict(set)
        
    def set_reference_length(self, length: int):
        """Set the reference length for relative position calculations."""
        self.reference_length = length
        
    def add_variant(self, variant_type: str, position: int, length: int, success: bool, affected_lineages: set):
        """Add a variant to tracking, handling SNPs and None types properly."""
        # Skip None mutations or mutations without proper data
        if variant_type is None or variant_type == "SKIPPED":
            return
            
        # Skip unsuccessful mutations
        if not success:
            return
            
        # Skip mutations without a position
        if position is None:
            return
        
        # Track variant for lineage analysis
        if affected_lineages:
            variant_key = (variant_type, position)
            self.variant_lineages[variant_key].update(affected_lineages)
            for lineage in affected_lineages:
                self.variants_by_lineage[lineage].add(variant_key)
        
        # Handle SNPs specially - they have no length but we want to track them
        if variant_type == "SNP":
            # Use size 1 for SNPs for visualization purposes
            self.variant_sizes[variant_type].append(1)
            self.variant_positions[variant_type].append((position, 1))
            self.max_position = max(self.max_position, position)
        elif length is not None and length > 0:
            # Handle other variants with actual lengths
            self.variant_sizes[variant_type].append(length)
            self.variant_positions[variant_type].append((position, length))
            self.max_position = max(self.max_position, position)
    
    def set_lineage_lengths(self, lineage_lengths: dict):
        self.lineage_lengths = lineage_lengths
        # If reference length not set, use the ancestral or first lineage length
        if not self.reference_length:
            if "Ancestral" in lineage_lengths:
                self.reference_length = lineage_lengths["Ancestral"]
            else:
                self.reference_length = max(lineage_lengths.values())

    def _get_relative_position(self, position: int) -> float:
        """Convert absolute position to relative position (0-100%)."""
        if self.reference_length and self.reference_length > 0:
            return (position / self.reference_length) * 100
        return position  # Fallback to absolute if no reference

    def save_variant_density_plot(self, output_path: str):
        """Create a density plot showing variant counts by RELATIVE position for each type."""
        MScompute(f"Creating variant density plot for sample {self.sample} chr {self.chromosome}")
        if not self.variant_positions:
            MSwarning("Warning: No variant data to plot")
            return
            
        plt.style.use('seaborn-v0_8-darkgrid')
        
        # Create subplots for each variant type
        variant_types = sorted(self.variant_positions.keys())
        n_types = len(variant_types)
        
        fig, axes = plt.subplots(n_types, 1, figsize=(14, 2.5 * n_types), sharex=True)
        if n_types == 1:
            axes = [axes]
        
        colors = {
            'SNP': '#e91e63',
            'INS': '#2ecc71',
            'DEL': '#e74c3c',
            'INV': '#3498db',
            'DUP': '#f39c12',
            'REPL': '#9b59b6'
        }
        
        for ax, var_type in zip(axes, variant_types):
            # Convert to relative positions
            positions = [self._get_relative_position(pos) for pos, _ in self.variant_positions[var_type]]
            
            # Create histogram with many bins for high resolution
            n_bins = min(100, len(positions) // 10)  # Use 100 bins for percentage scale
            if n_bins < 20:
                n_bins = 20
            
            ax.hist(positions, bins=n_bins, range=(0, 100),
                   color=colors.get(var_type, '#95a5a6'), 
                   alpha=0.7, edgecolor='none')
            
            ax.set_ylabel(f'{var_type}\nCount', fontsize=10)
            ax.set_xlim(0, 100)
            ax.grid(True, alpha=0.3)
            
            # Add count annotation
            ax.text(0.98, 0.95, f'n = {len(positions):,}', 
                   transform=ax.transAxes, ha='right', va='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        axes[-1].set_xlabel('Relative Genomic Position (%)', fontsize=12)
        fig.suptitle(f'Variant Density by Type - {self.sample} Chr{self.chromosome}',
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def save_size_distribution_plot(self, output_path: str):
        """Size distribution with SMOOTH DENSITY LINES instead of histograms."""
        MScompute(f"Creating variant size distribution for sample {self.sample} chr {self.chromosome}")
        
        # Filter out SNPs for the size distribution plot
        variant_sizes_no_snp = {k: v for k, v in self.variant_sizes.items() if k != 'SNP'}
        
        if not variant_sizes_no_snp:
            MSwarning("No non-SNP variant data to plot")
            return
            
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, ax = plt.subplots(figsize=(12, 8))

        colors = {
            'INS': '#2ecc71',
            'DEL': '#e74c3c',
            'INV': '#3498db',
            'DUP': '#f39c12',
            'REPL': '#9b59b6'
        }
        
        # Plot smooth density lines for each variant type
        for var_type in sorted(variant_sizes_no_snp.keys()):
            sizes = variant_sizes_no_snp[var_type]
            if len(sizes) > 1:  # Need at least 2 points for KDE
                # Create log-space for x-axis
                log_sizes = np.log10(sizes)
                
                try:
                    # Use kernel density estimation for smooth curve
                    kde = gaussian_kde(log_sizes, bw_method='scott')
                    
                    # Generate smooth x values in log space
                    x_min, x_max = np.min(log_sizes), np.max(log_sizes)
                    x_smooth = np.linspace(x_min - 0.5, x_max + 0.5, 500)
                    
                    # Calculate density
                    density = kde(x_smooth)
                    
                    # Convert back to linear space for plotting
                    x_linear = 10 ** x_smooth
                    
                    # Scale density by number of observations for better comparison
                    density_scaled = density * len(sizes)
                    
                    # Plot smooth line
                    ax.plot(x_linear, density_scaled,
                           label=f'{var_type} (n={len(sizes):,})',
                           color=colors.get(var_type, '#95a5a6'),
                           linewidth=2.5,
                           alpha=0.8)
                    
                    # Add shaded area under curve
                    ax.fill_between(x_linear, density_scaled, 
                                   alpha=0.2, 
                                   color=colors.get(var_type, '#95a5a6'))
                    
                except Exception as e:
                    MSwarning(f"Could not create smooth curve for {var_type}: {e}")
                    # Fallback to histogram if KDE fails
                    bins = np.logspace(np.log10(min(sizes)), np.log10(max(sizes)), 30)
                    ax.hist(sizes, bins=bins, 
                           label=f'{var_type} (n={len(sizes):,})',
                           color=colors.get(var_type, '#95a5a6'),
                           alpha=0.6,
                           edgecolor='black',
                           linewidth=0.5)
            elif len(sizes) == 1:
                # Single point - show as vertical line
                ax.axvline(x=sizes[0], 
                          label=f'{var_type} (n=1)',
                          color=colors.get(var_type, '#95a5a6'),
                          linewidth=2,
                          linestyle='--',
                          alpha=0.8)

        ax.set_xlabel('Variant Size (bp)', fontsize=12)
        ax.set_ylabel('Density (scaled by count)', fontsize=12)
        ax.set_xscale('log')
        ax.set_title(f'Variant Size Distribution - {self.sample} Chr{self.chromosome}',
                    fontsize=14, fontweight='bold')
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, which='both')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def save_cumulative_variants_plot(self, output_path: str):
        """Create a cumulative plot showing variant accumulation along RELATIVE positions."""
        MScompute(f"Creating cumulative variants plot for sample {self.sample} chr {self.chromosome}")
        
        if not self.variant_positions:
            MSwarning("No variant data to plot")
            return
        
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, ax = plt.subplots(figsize=(14, 8))
        
        colors = {
            'SNP': '#e91e63',
            'INS': '#2ecc71',
            'DEL': '#e74c3c',
            'INV': '#3498db',
            'DUP': '#f39c12',
            'REPL': '#9b59b6'
        }
        
        # Plot cumulative count for each variant type
        for var_type in sorted(self.variant_positions.keys()):
            # Convert to relative positions and sort
            relative_positions = sorted([self._get_relative_position(pos) 
                                        for pos, _ in self.variant_positions[var_type]])
            cumulative_counts = np.arange(1, len(relative_positions) + 1)
            
            ax.plot(relative_positions, cumulative_counts, 
                   color=colors.get(var_type, '#95a5a6'),
                   label=f'{var_type} (n={len(relative_positions):,})',
                   linewidth=2, alpha=0.8)
        
        ax.set_xlabel('Relative Genomic Position (%)', fontsize=12)
        ax.set_ylabel('Cumulative Variant Count', fontsize=12)
        ax.set_title(f'Cumulative Variant Distribution - {self.sample} Chr{self.chromosome}',
                    fontsize=14, fontweight='bold')
        ax.legend(loc='upper left', framealpha=0.9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 100)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    # Keep other methods unchanged
    def save_shared_variants_heatmap(self, output_path: str):
        """Create a heatmap showing shared variants between lineages."""
        MScompute(f"Creating shared variants heatmap for sample {self.sample} chr {self.chromosome}")
        
        if not self.variants_by_lineage:
            MSwarning("No variant data for lineage analysis")
            return
        
        # Get all lineages (excluding Ancestral if present)
        lineages = sorted([l for l in self.variants_by_lineage.keys() if l != "Ancestral"])
        
        if len(lineages) < 2:
            MSwarning("Need at least 2 lineages for shared variant analysis")
            return
        
        # Create matrix of shared variants
        n = len(lineages)
        shared_matrix = np.zeros((n, n))
        
        for i, lineage1 in enumerate(lineages):
            variants1 = self.variants_by_lineage[lineage1]
            for j, lineage2 in enumerate(lineages):
                if i == j:
                    # Diagonal: total variants in this lineage
                    shared_matrix[i, j] = len(variants1)
                else:
                    variants2 = self.variants_by_lineage[lineage2]
                    # Off-diagonal: shared variants
                    shared_matrix[i, j] = len(variants1.intersection(variants2))
        
        # Create heatmap
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(10, 8))
        
        im = ax.imshow(shared_matrix, cmap='YlOrRd', aspect='auto')
        
        # Set ticks
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))
        ax.set_xticklabels([f'L{l}' for l in lineages], rotation=45, ha='right')
        ax.set_yticklabels([f'L{l}' for l in lineages])
        
        # Add text annotations
        for i in range(n):
            for j in range(n):
                text = ax.text(j, i, f'{int(shared_matrix[i, j])}',
                             ha="center", va="center", 
                             color="white" if shared_matrix[i, j] > shared_matrix.max()/2 else "black",
                             fontsize=8)
        
        ax.set_title(f'Shared Variants Between Lineages - {self.sample} Chr{self.chromosome}',
                    fontsize=14, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Number of Variants', rotation=270, labelpad=20)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def save_variant_type_proportions_plot(self, output_path: str):
        """Create a pie chart showing proportions of each variant type."""
        MScompute(f"Creating variant type proportions plot for sample {self.sample} chr {self.chromosome}")
        
        if not self.variant_positions:
            MSwarning("Warning: No variant data to plot")
            return
        
        plt.style.use('default')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        colors = {
            'SNP': '#e91e63',
            'INS': '#2ecc71',
            'DEL': '#e74c3c',
            'INV': '#3498db',
            'DUP': '#f39c12',
            'REPL': '#9b59b6'
        }
        
        # Count variants by type
        variant_counts = {vtype: len(positions) for vtype, positions in self.variant_positions.items()}
        
        # Pie chart
        types = list(variant_counts.keys())
        counts = list(variant_counts.values())
        pie_colors = [colors.get(t, '#95a5a6') for t in types]
        
        wedges, texts, autotexts = ax1.pie(counts, labels=types, colors=pie_colors, 
                                           autopct='%1.1f%%', startangle=90)
        
        # Make percentage text bold
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_weight('bold')
        
        ax1.set_title('Variant Type Proportions', fontsize=12, fontweight='bold')
        
        # Bar chart for absolute numbers
        ax2.bar(types, counts, color=pie_colors, edgecolor='black', linewidth=1)
        ax2.set_ylabel('Number of Variants', fontsize=12)
        ax2.set_xlabel('Variant Type', fontsize=12)
        ax2.set_title('Variant Counts by Type', fontsize=12, fontweight='bold')
        ax2.grid(True, axis='y', alpha=0.3)
        
        # Add count labels on bars
        for i, (t, c) in enumerate(zip(types, counts)):
            ax2.text(i, c + max(counts)*0.01, f'{c:,}', 
                    ha='center', va='bottom', fontsize=10)
        
        fig.suptitle(f'Variant Type Analysis - {self.sample} Chr{self.chromosome}',
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

    def save_lineage_lengths_plot(self, output_path: str):
        """Original lineage lengths plot - kept as is."""
        MScompute(f"Creating lineage plot for sample {self.sample} chr {self.chromosome}")
        if not self.lineage_lengths:
            MSwarning("Warning: No lineage length data to plot")
            return
        
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, ax = plt.subplots(figsize=(10, max(6, len(self.lineage_lengths) * 0.3)))
        
        sorted_lineages = sorted(self.lineage_lengths.items(), key=lambda x: (isinstance(x[0], str), x[0]))
        lineages, lengths = zip(*sorted_lineages)
        
        lineage_labels = [f"Lineage {l}" if l != "Ancestral" else "Ancestral" for l in lineages]
        y_positions = np.arange(len(lineages))
        bars = ax.barh(y_positions, lengths, 
                       color='#3498db',
                       alpha=0.7,
                       edgecolor='black',
                       linewidth=1)

        for i, (bar, length) in enumerate(zip(bars, lengths)):
            ax.text(length + max(lengths) * 0.01, bar.get_y() + bar.get_height()/2,
                   f'{length:,} bp',
                   ha='left', va='center', fontsize=9)
        
        ax.set_yticks(y_positions)
        ax.set_yticklabels(lineage_labels)
        ax.set_xlabel('Sequence Length (bp)', fontsize=12)
        ax.set_title(f'Final Lineage Lengths - {self.sample} Chr{self.chromosome}',
                     fontsize=14, fontweight='bold')
        ax.grid(True, axis='x', alpha=0.3)
        ax.set_xlim(0, max(lengths) * 1.15)

        mean_length = np.mean(lengths)
        std_length = np.std(lengths)
        non_ancestral_lengths = [length for lineage, length in sorted_lineages if lineage != "Ancestral"]
        
        if non_ancestral_lengths:
            mean_non_ancestral = np.mean(non_ancestral_lengths)
            std_non_ancestral = np.std(non_ancestral_lengths)
            stats_text = (f"All lineages: μ={mean_length:,.0f} ± {std_length:,.0f} bp\n"
                          f"Non-ancestral: μ={mean_non_ancestral:,.0f} ± {std_non_ancestral:,.0f} bp")
        else:
            stats_text = f"Mean: {mean_length:,.0f} ± {std_length:,.0f} bp"
        
        ax.text(0.98, 0.02, stats_text,
                transform=ax.transAxes,
                ha='right', va='bottom',
                fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
class LintVisualizer:
    def __init__(self):
        self.total_nodes_before = 0
        self.total_nodes_after = 0
        self.removed_nodes = set()  # Names of removed nodes

    def record(self, before: int, after: int, removed: set):
        self.total_nodes_before = before
        self.total_nodes_after = after
        self.removed_nodes = removed

    def write_txt_report(self, output_path: str):
        try:
            num_removed = len(self.removed_nodes)
            percent_removed = (
                (num_removed / self.total_nodes_before) * 100
                if self.total_nodes_before else 0
            )

            with open(output_path, "a") as f:
                f.write("\n"+"-"*120+"\n")
                f.write("Lint Summary Report\n")
                f.write(f"Total nodes: {self.total_nodes_before}\n")
                f.write(f"Remaining nodes after linting: {self.total_nodes_after}\n")
                f.write(f"Orphan nodes removed by linting: {num_removed}\n")
                f.write(f"Percent removed: {percent_removed:.2f}%\n\n")
            
            MScompute(f"Saved lint stats and removed node list to {output_path}")
        except Exception as e:
            MSwarning(f"Could not write lint report: {e}")
