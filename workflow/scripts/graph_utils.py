"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

This scripts holds usefull function for the graph creation step 
"""

import random
from io_handler import MSerror, MSsuccess, MSwarning, MScompute
from datetime import datetime
import matplotlib.pyplot as plt # type: ignore
import matplotlib.patches as mpatches # type: ignore
import numpy as np
from collections import defaultdict

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
                
                length_str = str(mut['length']) if mut['length'] is not None else "N/A"
                status = 'SUCCESS' if mut['success'] else 'FAILED'
                
                f.write(f"{i:<12}\t{mut['tree_index']:<12}\t{mut['node_id']:<12}\t"
                       f"{mut['type']:<8}\t{mut['position']:<10}\t{length_str:<8}\t"
                       f"{lineages_str:<30}\t{status:<10}")
                
                # Add error message on the same line if failed
                if mut['error'] and not mut['success']:
                    f.write(f"\t# {mut['error']}")
                
                f.write("\n")
            
            f.write("-" * 120 + "\n")
            f.write(f"End of recap - Total mutations: {len(self.mutations)}\n")

class VariantSizeVisualizer:
    """Creates visualization graphs for variant sizes from mutation data."""
    
    def __init__(self, sample: str, chromosome: str):
        self.sample = sample
        self.chromosome = chromosome
        self.variant_sizes = defaultdict(list)
        self.variant_positions = defaultdict(list)
        self.lineage_lengths = {}
        self.max_position = 0
        
    def add_variant(self, variant_type: str, position: int, length: int, success: bool):
        if variant_type == "SNP" or length is None or not success:
            return
        self.variant_sizes[variant_type].append(length)
        self.variant_positions[variant_type].append((position, length))
        self.max_position = max(self.max_position, position)
    
    def set_lineage_lengths(self, lineage_lengths: dict):
        self.lineage_lengths = lineage_lengths

    def save_size_by_position_plot(self, output_path: str):
        MScompute("Creating variant size plot for sample {self.sample} chr {self.chromosome}")
        if not self.variant_positions:
            MSwarning("No variant data to plot (excluding SNPs)")
            return
            
        plt.style.use('seaborn-v0_8-darkgrid')
        fig, ax = plt.subplots(figsize=(14, 8))

        colors = {
            'INS': '#2ecc71',
            'DEL': '#e74c3c',
            'INV': '#3498db',
            'DUP': '#f39c12',
            'REPL': '#9b59b6'
        }
        
        for var_type in sorted(self.variant_positions.keys()):
            positions, sizes = zip(*self.variant_positions[var_type])
            ax.scatter(positions, sizes, 
                       c=colors.get(var_type, '#95a5a6'),
                       label=var_type, 
                       alpha=0.6, 
                       s=20,      
                       edgecolors='black',
                       linewidth=0.5)
        
        x_max = max(self.lineage_lengths.values()) if self.lineage_lengths else self.max_position * 1.05
        ax.set_xlim(0, x_max)

        # Ensure linear scale
        ax.set_yscale('linear')
        
        ax.set_xlabel('Genomic Position', fontsize=12)
        ax.set_ylabel('Variant Size (bp)', fontsize=12)
        ax.set_title(f'Variant Sizes by Position - {self.sample} Chr{self.chromosome}',
                     fontsize=14, fontweight='bold')
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, which='both')
        ax.grid(True, which='minor', alpha=0.1)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=1000, bbox_inches='tight')
        plt.close()

    
    def save_size_distribution_plot(self, output_path: str):
        MScompute(f"Creating variant size distribution for sample {self.sample} chr {self.chromosome}")
        if not self.variant_sizes:
            MSwarning("No variant data to plot (excluding SNPs)")
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
        
        all_sizes = []
        for sizes in self.variant_sizes.values():
            all_sizes.extend(sizes)
        all_sizes = [s for s in all_sizes if s > 0]
        if not all_sizes:
            MSwarning("No positive variant sizes available for plotting.")
            return

        min_size = min(all_sizes)
        max_size = max(all_sizes)
        bins = np.linspace(min_size, max_size, 50)

        sorted_types = sorted(self.variant_sizes.keys())
        hist_data = []
        colors_to_use = []

        for var_type in sorted_types:
            hist_data.append(self.variant_sizes[var_type])
            colors_to_use.append(colors.get(var_type, '#95a5a6'))
        
        ax.hist(hist_data, bins=bins, 
                label=sorted_types,
                color=colors_to_use,
                stacked=True,
                alpha=0.8,
                edgecolor='black',
                linewidth=0.5)

        ax.set_xlabel('Variant Size (bp)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'Variant Size Distribution - {self.sample} Chr{self.chromosome}',
                     fontsize=14, fontweight='bold')
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3, which='both')

        stats_text = []
        for var_type in sorted_types:
            sizes = self.variant_sizes[var_type]
            if sizes:
                stats_text.append(
                    f"{var_type}: n={len(sizes)}, "
                    f"median={np.median(sizes):.0f}bp, "
                    f"range={min(sizes)}-{max(sizes)}bp"
                )
        
        ax.text(0.02, 0.98, '\n'.join(stats_text), 
                transform=ax.transAxes,
                verticalalignment='top', 
                fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

        plt.tight_layout()
        plt.savefig(output_path, dpi=1000, bbox_inches='tight')
        plt.close()


    def save_lineage_lengths_plot(self, output_path: str):
        MScompute(f"Creating lineage plot for sample {self.sample} chr {self.chromosome}")
        if not self.lineage_lengths:
            MSwarning("No lineage length data to plot")
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
        plt.savefig(output_path, dpi=1000, bbox_inches='tight')
        plt.close()
