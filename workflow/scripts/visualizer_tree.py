#!/usr/bin/env python3
import os
import tskit
import random
import matplotlib.pyplot as plt
import argparse
import sys
from io_handler import MSerror, MSsuccess, MScompute


class ProgressTracker:
    """Simple progress tracker for visualization tasks"""
    def __init__(self, total_tasks):
        self.total_tasks = total_tasks
        self.completed_tasks = 0
        self.current_subtask = 0
        self.total_subtasks = 0
    
    def start_task(self, task_name, subtasks=1):
        self.current_subtask = 0
        self.total_subtasks = subtasks
        MScompute(f"{task_name}... [{self.get_overall_progress():.0f}%]")
    
    def update_subtask(self):
        self.current_subtask += 1
        if self.current_subtask == self.total_subtasks:
            self.completed_tasks += 1
    
    def get_overall_progress(self):
        return (self.completed_tasks / self.total_tasks) * 100


class TreeVisualizer:
    def __init__(self, mutated_ts, ancestry_ts, output_dir="visualizations", edge_width=3):
        self.mutated_ts = mutated_ts
        self.ancestry_ts = ancestry_ts
        self.output_dir = output_dir
        self.edge_width = edge_width
        os.makedirs(output_dir, exist_ok=True)
        
        # Base CSS for thicker edges
        self.base_style = f".edge {{stroke-width: {self.edge_width}px}}"
        
        # Initialize progress tracker
        self.progress = ProgressTracker(total_tasks=4)  # Adjust based on tasks
    
    def save_global_trees_svg(self):
        """Save global tree visualizations as SVG"""
        self.progress.start_task("Saving global tree SVGs", 2)
        
        # Mutated tree sequence
        with open(os.path.join(self.output_dir, "mutated_global.svg"), "w") as f:
            f.write(self.mutated_ts.draw_svg(time_scale="rank", style=self.base_style))
        self.progress.update_subtask()
        
        # Ancestry tree sequence
        with open(os.path.join(self.output_dir, "ancestry_global.svg"), "w") as f:
            f.write(self.ancestry_ts.draw_svg(time_scale="rank", style=self.base_style))
        self.progress.update_subtask()
    
    def save_individual_trees(self, chromosome_name="chr1"):
        """Save ALL individual tree visualizations"""
        total_trees = self.mutated_ts.num_trees + self.ancestry_ts.num_trees
        self.progress.start_task(f"Saving {total_trees} individual trees", total_trees)
        
        # Create subdirectories
        mutated_dir = os.path.join(self.output_dir, "individual_trees", "mutated")
        ancestry_dir = os.path.join(self.output_dir, "individual_trees", "ancestry")
        os.makedirs(mutated_dir, exist_ok=True)
        os.makedirs(ancestry_dir, exist_ok=True)
        
        # Save ALL mutated trees
        for i, tree in enumerate(self.mutated_ts.trees()):
            filename = f"chr_{chromosome_name}_mutated_tree_{i}.svg"
            path = os.path.join(mutated_dir, filename)
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
            self.progress.update_subtask()
        
        # Save ALL ancestry trees
        for i, tree in enumerate(self.ancestry_ts.trees()):
            filename = f"chr_{chromosome_name}_ancestry_tree_{i}.svg"
            path = os.path.join(ancestry_dir, filename)
            
            with open(path, "w") as f:
                f.write(tree.draw_svg(time_scale="rank", style=self.base_style))
            self.progress.update_subtask()
    
    def save_highlighted_trees_same_branch(self, num_occasions=5):
        """Save highlighted trees where the same branch is highlighted across all trees in each occasion"""
        total_operations = num_occasions * self.ancestry_ts.num_trees
        self.progress.start_task(f"Creating highlighted trees ({num_occasions} occasions)", total_operations)
        
        highlight_base_dir = os.path.join(self.output_dir, "highlighted_trees")
        os.makedirs(highlight_base_dir, exist_ok=True)
        
        total_trees = self.ancestry_ts.num_trees
        highlight_width = self.edge_width + 3
        
        for occasion in range(num_occasions):
            # Create subdirectory for this occasion
            occasion_dir = os.path.join(highlight_base_dir, f"occasion_{occasion + 1}")
            os.makedirs(occasion_dir, exist_ok=True)
            
            # Pick a random node to highlight consistently across all trees
            first_tree = self.ancestry_ts.first()
            nodes = list(first_tree.nodes())
            
            if len(nodes) > 1:
                internal_nodes = [n for n in nodes if first_tree.num_children(n) > 0]
                if internal_nodes:
                    target_node = random.choice(internal_nodes)
                else:
                    target_node = random.choice(nodes[:-1])
                
                # Choose a consistent style for this occasion
                styles = [
                    ("cyan", "5,5", "Cyan dashed"),
                    ("red", "none", "Red solid"),
                    ("orange", "10,5", "Orange dashed"),
                    ("purple", "2,2", "Purple dotted"),
                    ("green", "8,3", "Green dashed")
                ]
                
                color, dash, description = styles[occasion % len(styles)]
                
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
                        
                        self.progress.update_subtask()
                                
                    except Exception as e:
                        self.progress.update_subtask()
                        continue
                
                # Save occasion info
                info_content = f"Occasion {occasion + 1}: Highlighting node {target_node}\n"
                info_content += f"Style: {description}\n"
                info_content += f"Total trees processed: {total_trees}\n"
                info_content += f"Trees with node highlighted: {trees_saved}\n"
                info_content += f"Trees without node: {total_trees - trees_saved}\n"
                
                with open(os.path.join(occasion_dir, "occasion_info.txt"), "w") as f:
                    f.write(info_content)
    
    def save_additional_visualizations(self):
        """Additional visualizations inspired by tskit tutorial"""
        self.progress.start_task("Creating additional visualizations", 3)
        
        extra_dir = os.path.join(self.output_dir, "extra_visualizations")
        os.makedirs(extra_dir, exist_ok=True)
        
        # 1. Time scale visualizations (first tree only for efficiency)
        time_scales = ["time", "rank"]
        for scale in time_scales:
            try:
                first_tree = self.ancestry_ts.first()
                filename = f"first_tree_{scale}_scale.svg"
                path = os.path.join(extra_dir, filename)
                
                with open(path, "w") as f:
                    f.write(first_tree.draw_svg(time_scale=scale, size=(800, 600), style=self.base_style))
            except Exception as e:
                pass
        self.progress.update_subtask()
        
        # 2. Tree with custom node labels (first tree only)
        try:
            first_tree = self.ancestry_ts.first()
            node_labels = {node.id: f"n{node.id}" for node in first_tree.nodes()}
            with open(os.path.join(extra_dir, "first_tree_with_node_labels.svg"), "w") as f:
                f.write(first_tree.draw_svg(node_labels=node_labels, size=(800, 600), style=self.base_style))
        except Exception as e:
            pass
        self.progress.update_subtask()
        
        # 3. Tree sequence overview plot
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
            plt.savefig(os.path.join(extra_dir, "tree_sequence_overview.png"), dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            pass
        self.progress.update_subtask()


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
    
    parser.add_argument("mutated_ts", help="Path to the mutated tree sequence file (.trees format)")
    parser.add_argument("ancestry_ts", help="Path to the ancestry tree sequence file (.trees format)")
    parser.add_argument("-o", "--output-dir", default="poster_visualizations", help="Output directory for visualizations (default: poster_visualizations)")
    parser.add_argument("-c", "--chromosome", default="chr1", help="Chromosome name for file naming (default: chr1)")
    parser.add_argument("--occasions", type=int, default=5, help="Number of highlighting occasions (default: 5)")
    parser.add_argument("--edge-width", type=int, default=3, help="Base edge width in pixels (default: 3)")
    parser.add_argument("--skip-individual", action="store_true", help="Skip saving individual tree SVGs (useful for large datasets)")
    parser.add_argument("--skip-extra", action="store_true", help="Skip additional visualizations")
    
    return parser.parse_args()


def main():
    """Main function to run all visualizations"""
    args = parse_arguments()
    
    # Validate input files exist
    for filepath in [args.mutated_ts, args.ancestry_ts]:
        if not os.path.exists(filepath):
            MSerror(f"Error: File does not exist: {filepath}")
            sys.exit(1)
    
    MScompute("Loading tree sequences...")
    mutated_ts, ancestry_ts = load_tree_sequences(args.mutated_ts, args.ancestry_ts)
    
    # Create visualizer with custom edge width
    visualizer = TreeVisualizer(mutated_ts, ancestry_ts, args.output_dir, edge_width=args.edge_width)
    
    # Generate visualizations based on arguments
    visualizer.save_global_trees_svg()
    
    if not args.skip_individual:
        visualizer.save_individual_trees(chromosome_name=args.chromosome)
    
    visualizer.save_highlighted_trees_same_branch(num_occasions=args.occasions)
    
    if not args.skip_extra:
        visualizer.save_additional_visualizations()
    
    MSsuccess("All visualizations completed! [100%]")


if __name__ == "__main__":
    main()