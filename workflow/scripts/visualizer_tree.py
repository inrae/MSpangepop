#!/usr/bin/env python3

"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage: Produce simplified visualizations for the tree sequence 
"""

import os
import tskit
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
import matplotlib.pyplot as plt
import argparse
import sys
from io_handler import MSerror, MSsuccess, MScompute


class TreeVisualizer:
    def __init__(self, mutated_ts, ancestry_ts, output_dir="visualizations", edge_width=3):
        self.mutated_ts = mutated_ts
        self.ancestry_ts = ancestry_ts
        self.output_dir = output_dir
        self.edge_width = edge_width
        os.makedirs(output_dir, exist_ok=True)
        
        # Base CSS for thicker edges
        self.base_style = f".edge {{stroke-width: {self.edge_width}px}}"
    
    def save_global_trees_svg(self):
        """Save global tree visualizations as SVG"""
        MScompute("Saving global tree SVGs...")
        
        # Mutated tree sequence
        with open(os.path.join(self.output_dir, "mutated_global.svg"), "w") as f:
            f.write(self.mutated_ts.draw_svg(style=self.base_style))
        
        # Ancestry tree sequence
        with open(os.path.join(self.output_dir, "ancestry_global.svg"), "w") as f:
            f.write(self.ancestry_ts.draw_svg(style=self.base_style))
    
    def save_individual_trees(self, chromosome_name="chr1", max_trees=100):
        """Save first 100 individual tree visualizations"""
        MScompute(f"Saving individual trees (up to {max_trees})...")
        
        # Create subdirectories
        mutated_dir = os.path.join(self.output_dir, "individual_trees", "mutated")
        ancestry_dir = os.path.join(self.output_dir, "individual_trees", "ancestry")
        os.makedirs(mutated_dir, exist_ok=True)
        os.makedirs(ancestry_dir, exist_ok=True)
        
        # Save mutated trees
        for i, tree in enumerate(self.mutated_ts.trees()):
            if i >= max_trees:
                break
                
            filename = f"chr_{chromosome_name}_mutated_tree_{i}.svg"
            path = os.path.join(mutated_dir, filename)
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))

        # Save ancestry trees
        for i, tree in enumerate(self.ancestry_ts.trees()):
            if i >= max_trees:
                break
                
            filename = f"chr_{chromosome_name}_ancestry_tree_{i}.svg"
            path = os.path.join(ancestry_dir, filename)
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
    
    def save_tree_sequence_overview(self):
        """Tree sequence overview plot"""
        MScompute("Creating tree sequence overview...")
        
        try:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            positions = []
            tree_indices = []
            
            for i, tree in enumerate(self.ancestry_ts.trees()):
                positions.extend([tree.interval.left, tree.interval.right])
                tree_indices.extend([i, i])
            
            ax.plot(positions, tree_indices, 'b-', linewidth=3)
            ax.set_xlabel('Genomic Position', fontsize=12)
            ax.set_ylabel('Tree Index', fontsize=12)
            ax.set_title('Tree Sequence Overview', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "tree_sequence_overview.png"), 
                       dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            MSerror(f"Error creating overview plot: {e}")


def load_tree_sequences(mutated_path, ancestry_path):
    """Load both tree sequences from specified file paths"""
    try:
        mutated_ts = tskit.load(mutated_path)
        ancestry_ts = tskit.load(ancestry_path)
        return mutated_ts, ancestry_ts
    except Exception as e:
        MSerror(f"Error loading tree sequences: {e}")
        sys.exit(1)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate simplified tree sequence visualizations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tree_visualizer.py mutated.trees ancestry.trees
  python tree_visualizer.py mutated.trees ancestry.trees -o my_visualizations
  python tree_visualizer.py mutated.trees ancestry.trees -c chr2 --edge-width 5
        """
    )
    
    parser.add_argument("mutated_ts", help="Path to the mutated tree sequence file (.trees format)")
    parser.add_argument("ancestry_ts", help="Path to the ancestry tree sequence file (.trees format)")
    parser.add_argument("-o", "--output-dir", default="visualizations", 
                       help="Output directory for visualizations (default: visualizations)")
    parser.add_argument("-c", "--chromosome", default="chr1", 
                       help="Chromosome name for file naming (default: chr1)")
    parser.add_argument("--edge-width", type=int, default=3, 
                       help="Base edge width in pixels (default: 3)")
    
    return parser.parse_args()


def main():
    """Main function to run all visualizations"""
    args = parse_arguments()
    
    # Validate input files exist
    for filepath in [args.mutated_ts, args.ancestry_ts]:
        if not os.path.exists(filepath):
            MSerror(f"Error: File does not exist: {filepath}")
            sys.exit(1)
    
    MScompute("Starting to compute tree visualizations")
    mutated_ts, ancestry_ts = load_tree_sequences(args.mutated_ts, args.ancestry_ts)
    
    # Create visualizer
    visualizer = TreeVisualizer(mutated_ts, ancestry_ts, args.output_dir, edge_width=args.edge_width)
    
    # Generate visualizations
    visualizer.save_global_trees_svg()
    visualizer.save_individual_trees(chromosome_name=args.chromosome, max_trees=100)
    visualizer.save_tree_sequence_overview()
    
    MSsuccess("All tree visualizations completed!")


if __name__ == "__main__":
    main()