import os
import tskit
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import sys
from io_handler import MSerror, MSsuccess, MScompute, MSwarning

try:
    import tskit_arg_visualizer
    HAS_D3ARG = True
except ImportError:
    HAS_D3ARG = False
    MSwarning("Warning: tskit_arg_visualizer not available. Skipping D3ARG visualizations.")

class TreeVisualizer:
    def __init__(self, mutated_ts, ancestry_ts, output_dir="visualizations", edge_width=3):
        self.mutated_ts = mutated_ts
        self.ancestry_ts = ancestry_ts
        self.output_dir = output_dir
        self.edge_width = edge_width  # Default edge width
        os.makedirs(output_dir, exist_ok=True)
        
        # Base CSS for thicker edges
        self.base_style = f".edge {{stroke-width: {self.edge_width}px}}"
    
    def save_global_trees_svg(self):
        """Save global tree visualizations as SVG"""
        MScompute("Saving global tree SVGs...")
        
        # Mutated tree sequence
        with open(os.path.join(self.output_dir, "mutated_global.svg"), "w") as f:
            f.write(self.mutated_ts.draw_svg(time_scale="rank", style=self.base_style))
        
        # Ancestry tree sequence
        with open(os.path.join(self.output_dir, "ancestry_global.svg"), "w") as f:
            f.write(self.ancestry_ts.draw_svg(time_scale="rank", style=self.base_style))
        
    
    def save_individual_trees(self, chromosome_name="chr1"):
        """Save ALL individual tree visualizations"""
        MScompute("Saving ALL individual trees...")
        
        # Create subdirectories
        mutated_dir = os.path.join(self.output_dir, "individual_trees", "mutated")
        ancestry_dir = os.path.join(self.output_dir, "individual_trees", "ancestry")
        os.makedirs(mutated_dir, exist_ok=True)
        os.makedirs(ancestry_dir, exist_ok=True)
        
        # Save ALL mutated trees
        total_mutated_trees = self.mutated_ts.num_trees
        MScompute(f"Saving {total_mutated_trees} mutated trees...")
        
        for i, tree in enumerate(self.mutated_ts.trees()):
            filename = f"chr_{chromosome_name}_mutated_tree_{i}.svg"
            path = os.path.join(mutated_dir, filename)
            
            if i % 100 == 0 or i == total_mutated_trees - 1:
                MScompute(f"Saving mutated tree {i+1}/{total_mutated_trees}")
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
        
        # Save ALL ancestry trees
        total_ancestry_trees = self.ancestry_ts.num_trees
        MScompute(f"Saving {total_ancestry_trees} ancestry trees...")
        
        for i, tree in enumerate(self.ancestry_ts.trees()):
            filename = f"chr_{chromosome_name}_ancestry_tree_{i}.svg"
            path = os.path.join(ancestry_dir, filename)
            
            if i % 100 == 0 or i == total_ancestry_trees - 1:
                MScompute(f"Saving ancestry tree {i+1}/{total_ancestry_trees}")
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
        
    
    def save_highlighted_trees_same_branch(self, num_occasions=5):
        """Save highlighted trees where the same branch is highlighted across all trees in each occasion"""
        MScompute(f"Creating highlighted tree visualizations with same branch per occasion...")
        
        highlight_base_dir = os.path.join(self.output_dir, "highlighted_trees")
        os.makedirs(highlight_base_dir, exist_ok=True)
        
        total_trees = self.ancestry_ts.num_trees
        
        # Increased highlight width for better visibility
        highlight_width = self.edge_width + 3
        
        for occasion in range(num_occasions):
            MScompute(f"Creating highlighted trees for occasion {occasion + 1}/{num_occasions}")
            
            # Create subdirectory for this occasion
            occasion_dir = os.path.join(highlight_base_dir, f"occasion_{occasion + 1}")
            os.makedirs(occasion_dir, exist_ok=True)
            
            # Pick a random node to highlight consistently across all trees
            # We'll use the first tree to determine available nodes
            first_tree = self.ancestry_ts.first()
            nodes = list(first_tree.nodes())
            
            if len(nodes) > 1:
                # Try to get a node that exists in most trees
                internal_nodes = [n for n in nodes if first_tree.num_children(n) > 0]
                if internal_nodes:
                    target_node = random.choice(internal_nodes)
                else:
                    target_node = random.choice(nodes[:-1])  # Exclude root
                
                # Choose a consistent style for this occasion
                styles = [
                    ("cyan", "5,5", "Cyan dashed"),
                    ("red", "none", "Red solid"),
                    ("orange", "10,5", "Orange dashed"),
                    ("purple", "2,2", "Purple dotted"),
                    ("green", "8,3", "Green dashed")
                ]
                
                color, dash, description = styles[occasion % len(styles)]
                
                # Save info about this occasion
                info_content = f"Occasion {occasion + 1}: Highlighting node {target_node}\n"
                info_content += f"Style: {description}\n"
                info_content += f"Total trees processed: 0\n"  # Will be updated
                
                trees_saved = 0
                
                # Go through all trees and highlight the same node if it exists
                for tree_idx, tree in enumerate(self.ancestry_ts.trees()):
                    try:
                        # Check if the target node exists in this tree
                        if target_node in tree.nodes():
                            parent = tree.parent(target_node)
                            
                            if parent != tskit.NULL:
                                # Create CSS for highlighting with base style
                                if dash == "none":
                                    css_string = f"{self.base_style} .n{target_node} > .edge {{stroke: {color}; stroke-width: {highlight_width}px}}"
                                else:
                                    css_string = f"{self.base_style} .n{target_node} > .edge {{stroke: {color}; stroke-width: {highlight_width}px; stroke-dasharray: {dash}}}"
                                
                                # Save highlighted tree
                                filename = f"ancestry_tree_{tree_idx}_node_{target_node}_highlighted.svg"
                                path = os.path.join(occasion_dir, filename)
                                
                                with open(path, "w") as f:
                                    f.write(tree.draw_svg(time_scale="rank", style=css_string))
                                
                                trees_saved += 1
                        else:
                            # Node doesn't exist in this tree, save without highlighting but with base style
                            filename = f"ancestry_tree_{tree_idx}_node_{target_node}_not_present.svg"
                            path = os.path.join(occasion_dir, filename)
                            
                            with open(path, "w") as f:
                                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
                                
                    except Exception as e:
                        MSerror(f"Failed to process tree {tree_idx} for occasion {occasion + 1}: {e}")
                        continue
                
                # Update and save occasion info
                info_content = info_content.replace("Total trees processed: 0", f"Total trees processed: {total_trees}")
                info_content += f"Trees with node highlighted: {trees_saved}\n"
                info_content += f"Trees without node: {total_trees - trees_saved}\n"
                
                with open(os.path.join(occasion_dir, "occasion_info.txt"), "w") as f:
                    f.write(info_content)
                
                MScompute(f"Occasion {occasion + 1}: Highlighted node {target_node} in {trees_saved}/{total_trees} trees")
    
    def save_additional_visualizations(self):
        """Additional visualizations inspired by tskit tutorial"""
        MScompute("Creating additional visualizations...")
        
        extra_dir = os.path.join(self.output_dir, "extra_visualizations")
        os.makedirs(extra_dir, exist_ok=True)
        
        # 1. All trees with different time scales (time and rank for all trees)
        time_scales = ["time", "rank"]
        
        for scale in time_scales:
            MScompute(f"Creating {scale} scale visualizations for all trees...")
            
            scale_dir = os.path.join(extra_dir, f"all_trees_{scale}_scale")
            os.makedirs(scale_dir, exist_ok=True)
            
            try:
                # Save all ancestry trees with this time scale
                for i, tree in enumerate(self.ancestry_ts.trees()):
                    filename = f"ancestry_tree_{i}_{scale}_scale.svg"
                    path = os.path.join(scale_dir, filename)
                    
                    if i % 100 == 0 or i == self.ancestry_ts.num_trees - 1:
                        MScompute(f"Saving {scale} scale tree {i+1}/{self.ancestry_ts.num_trees}")
                    
                    with open(path, "w") as f:
                        f.write(tree.draw_svg(time_scale=scale, size=(800, 600), style=self.base_style))
                        
            except Exception as e:
                MSerror(f"Failed to create {scale} scale visualization: {e}")
        
        # 2. Tree with custom node labels (first tree only)
        try:
            first_tree = self.ancestry_ts.first()
            # Add node labels showing node IDs
            node_labels = {node.id: f"n{node.id}" for node in first_tree.nodes()}
            with open(os.path.join(extra_dir, "first_tree_with_node_labels.svg"), "w") as f:
                f.write(first_tree.draw_svg(node_labels=node_labels, size=(800, 600), style=self.base_style))
        except Exception as e:
            MSerror(f"Failed to create node labels visualization: {e}")
        
        # 3. Mutations visualization (if available)
        if self.mutated_ts.num_mutations > 0:
            try:
                first_tree_mut = self.mutated_ts.first()
                with open(os.path.join(extra_dir, "first_tree_with_mutations.svg"), "w") as f:
                    f.write(first_tree_mut.draw_svg(size=(800, 600), style=self.base_style))
            except Exception as e:
                MSerror(f"Failed to create mutations visualization: {e}")
        
        # 4. Tree sequence overview plot
        try:
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Plot tree spans with thicker lines
            positions = []
            tree_indices = []
            
            for i, tree in enumerate(self.ancestry_ts.trees()):
                positions.extend([tree.interval.left, tree.interval.right])
                tree_indices.extend([i, i])
            
            ax.plot(positions, tree_indices, 'b-', linewidth=3)  # Increased line width
            ax.set_xlabel('Genomic Position', fontsize=12)
            ax.set_ylabel('Tree Index', fontsize=12)
            ax.set_title('Tree Sequence Overview', fontsize=14)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(os.path.join(extra_dir, "tree_sequence_overview.png"), dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            MSerror(f"Failed to create overview plot: {e}")
        

def load_tree_sequences(mutated_path, ancestry_path):
    """Load both tree sequences from specified file paths"""
    MScompute(f"Loading mutated tree sequence")
    MScompute(f"Loading ancestry tree sequence")
    
    try:
        mutated_ts = tskit.load(mutated_path)
        ancestry_ts = tskit.load(ancestry_path)
        
        MScompute(f"Successfully loaded:")
        MScompute(f"  - Mutated TS: {mutated_ts.num_trees} trees, {mutated_ts.num_samples} samples")
        MScompute(f"  - Ancestry TS: {ancestry_ts.num_trees} trees, {ancestry_ts.num_samples} samples")
        
        return mutated_ts, ancestry_ts
        
    except Exception as e:
        MSerror(f"Error loading tree sequences: {e}")
        sys.exit(1)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate tree sequence visualizations for poster",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tree_visualizer.py mutated.trees ancestry.trees
  python tree_visualizer.py mutated.trees ancestry.trees -o my_visualizations
  python tree_visualizer.py mutated.trees ancestry.trees -c chr2 --occasions 3
  python tree_visualizer.py mutated.trees ancestry.trees --edge-width 5
        """
    )
    
    parser.add_argument(
        "mutated_ts",
        help="Path to the mutated tree sequence file (.trees format)"
    )
    
    parser.add_argument(
        "ancestry_ts", 
        help="Path to the ancestry tree sequence file (.trees format)"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        default="poster_visualizations",
        help="Output directory for visualizations (default: poster_visualizations)"
    )
    
    parser.add_argument(
        "-c", "--chromosome",
        default="chr1",
        help="Chromosome name for file naming (default: chr1)"
    )
    
    parser.add_argument(
        "--occasions",
        type=int,
        default=5,
        help="Number of highlighting occasions (default: 5)"
    )
    
    parser.add_argument(
        "--edge-width",
        type=int,
        default=3,
        help="Base edge width in pixels (default: 3)"
    )
    
    parser.add_argument(
        "--skip-individual",
        action="store_true",
        help="Skip saving individual tree SVGs (useful for large datasets)"
    )
    
    parser.add_argument(
        "--skip-extra",
        action="store_true",
        help="Skip additional visualizations"
    )
    
    return parser.parse_args()

def main():
    """Main function to run all visualizations"""
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Validate input files exist
    for filepath in [args.mutated_ts, args.ancestry_ts]:
        if not os.path.exists(filepath):
            MSerror(f"Error: File does not exist: {filepath}")
            sys.exit(1)
    
    # Load tree sequences
    mutated_ts, ancestry_ts = load_tree_sequences(args.mutated_ts, args.ancestry_ts)
    
    # Create visualizer with custom edge width
    visualizer = TreeVisualizer(mutated_ts, ancestry_ts, args.output_dir, edge_width=args.edge_width)
    
    MScompute(f"\nGenerating visualizations in: {args.output_dir}")
    MScompute(f"Using edge width: {args.edge_width}px")
    
    # Generate visualizations based on arguments
    visualizer.save_global_trees_svg()
    
    if not args.skip_individual:
        visualizer.save_individual_trees(chromosome_name=args.chromosome)
    else:
        MScompute("Skipping individual tree generation (--skip-individual specified)")
    
    visualizer.save_highlighted_trees_same_branch(num_occasions=args.occasions)
    
    if not args.skip_extra:
        visualizer.save_additional_visualizations()
    
    

if __name__ == "__main__":
    main()
    MSsuccess("Done !")