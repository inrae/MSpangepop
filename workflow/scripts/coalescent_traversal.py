import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from io_handler import MSpangepopDataHandler, MScompute, MSsuccess, MSerror

class TreeNode:
    """Represents a node in the binary tree."""
    
    def __init__(self, node_id, time):
        self.id = node_id  # Node ID
        self.time = time  # Time attribute for ordering nodes
        self.left = None
        self.right = None
        self.mutation_count = 0  # Store the count of mutations
        self.leaves = set()  # Set of all descendant leaf nodes, including itself if it is a leaf

class BinaryTree:
    """Handles tree construction and processing."""

    def __init__(self, tree_data):
        self.nodes = {}  # Dictionary to store TreeNode objects by ID
        self.root = None
        self.tree_index = tree_data.get("tree_index", "unknown")  # Ensure default value if missing
        self.interval = tree_data.get("interval", None)  # Retrieve interval if present
        self.build_tree(tree_data)  # Construct tree from input data

    def build_tree(self, tree_data):
        """Builds the binary tree from input JSON data with error handling."""
        try:
            # Step 1: Validate input structure
            if "nodes" not in tree_data or "edges" not in tree_data or "mutations" not in tree_data:
                raise MSerror(f"Malformed tree data for tree index {self.tree_index}. Missing required fields.")

            if self.interval is None:
                raise MSerror(f"Missing interval for tree index {self.tree_index}.")

            for node in tree_data["nodes"]:
                self.nodes[node["id"]] = TreeNode(node["id"], node["time"])

            # Step 2: Link nodes using parent-child relationships
            child_counts = defaultdict(int)  # Track how many children each parent has
            for edge in tree_data["edges"]:
                parent_id, child_id = edge["parent"], edge["child"]
                if parent_id not in self.nodes or child_id not in self.nodes:
                    raise MSerror(f"Tree structure inconsistency: Parent {parent_id} or Child {child_id} missing.")

                parent = self.nodes[parent_id]
                child = self.nodes[child_id]

                if child_counts[parent_id] == 0:
                    parent.left = child
                else:
                    parent.right = child
                child_counts[parent_id] += 1

            # Step 3: Assign mutation counts to nodes
            for mutation in tree_data["mutations"]:
                node_id = mutation["node"]
                if node_id not in self.nodes:
                    raise MSerror(f"Mutation assigned to non-existent node {node_id}.")
                self.nodes[node_id].mutation_count += 1  # Increment count instead of storing times

            # Step 4: Identify the root node
            all_children = {edge["child"] for edge in tree_data["edges"]}
            root_candidates = set(self.nodes.keys()) - all_children
            if root_candidates:
                self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]
            else:
                raise MSerror(f"Tree index {self.tree_index} has no identifiable root node.")

            # Step 5: Compute affected leaf sets
            self.compute_leaf_sets()

        except KeyError as e:
            raise MSerror(f"Key error during tree building for tree {self.tree_index}: {e}")
        except Exception as e:
            raise MSerror(f"Unexpected error while processing tree {self.tree_index}: {e}")

    def compute_leaf_sets(self):
        """Computes and stores the affected leaf nodes for each node in the tree."""
        
        def dfs(node):
            """Recursive DFS to populate affected leaf sets."""
            if not node:
                return set()
            
            # If a node is a leaf, it should contain only itself in the affected set
            if node.left is None and node.right is None:
                node.leaves = {node.id}
                return node.leaves

            # Recursively collect leaves
            left_leaves = dfs(node.left)
            right_leaves = dfs(node.right)

            # Store all descendant leaves
            node.leaves = left_leaves | right_leaves
            return node.leaves

        if not self.root:
            raise MSerror(f"Cannot compute leaf sets: No root node found for tree {self.tree_index}.")

        dfs(self.root)

    def get_mutated_nodes(self):
        """Returns only mutated nodes while keeping the tree structure."""
        mutated_nodes = []

        def traverse(node):
            """Recursively collect only mutated nodes."""
            if node:
                if node.mutation_count > 0:  # Only include nodes with mutations
                    mutated_nodes.append({
                        "node": node.id,
                        "mutation_count": node.mutation_count,  # Save the count of mutations
                        "affected_leaves": list(node.leaves)  # Keep affected leaves for reference
                    })
                traverse(node.left)
                traverse(node.right)

        traverse(self.root)
        return {
            "tree_index": self.tree_index,
            "interval": self.interval,  # Include interval in output
            "mutated_nodes": mutated_nodes
        }

def process_tree(tree_data):
    """
    Processes a single tree and returns mutated nodes with mutation count and interval.
    """
    try:
        tree = BinaryTree(tree_data)
        return tree.get_mutated_nodes()
    except MSerror as e:
        raise MSerror(f"Error while processing a tree: {e}")


def main(json_file, chromosome, output_file, threads, readable_json):
    """
    Main function to process tree data in parallel and save the results to a JSON file.
    
    Parameters:
        json_file (str): Path to the input JSON file containing tree data.
        chromosome (str): Chromosome identifier.
        output_file (str): Path to save the serialized traversal output.
        threads (int): Number of parallel threads to use.
    """
    
    try:
        # Load input data
        json_data = MSpangepopDataHandler.read_json(json_file)
        if not isinstance(json_data, list):
            raise MSerror(f"Invalid JSON format in {json_file}: Expected a list of tree objects.")

        MScompute(f"Processing mutated nodes of chromosome {chromosome} with {threads} threads")

        # Process trees in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=threads) as executor:
            results = list(filter(None, executor.map(process_tree, json_data)))  # Remove None values

        # Save output JSON, which now includes interval data
        if results:
            MSpangepopDataHandler.save_json(results, output_file, readable_json)
            MSsuccess(f"Serialized mutated nodes for chromosome {chromosome}")
        else:
            raise MSerror(f"No valid trees processed for chromosome {chromosome}. Output file not saved.")

    except FileNotFoundError:
        raise MSerror(f"Input file not found: {json_file}")
    except Exception as e:
        raise MSerror(f"Critical error during execution: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serialize mutated nodes into JSON with mutation counts and intervals.")
    parser.add_argument("--json", required=True, help="Path to the input JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome number")
    parser.add_argument("--output_file", required=True, help="Path to the output JSON file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False],
                        help="Save JSON in a human-readable format (True/False, default: False).")
    args = parser.parse_args()

    # Execute the main function
    try:
        main(args.json, args.chromosome, args.output_file, args.threads, args.readable_json)
    except MSerror as e:
        raise MSerror(f"Fatal error: {e}")
