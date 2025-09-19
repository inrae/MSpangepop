"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Construct a binary tree form a json file and then perform the perorder traversal on it.

--json Path to the JSON file containing the tree
--output_file Path to the output JSON 
--chromosome Chromosome name
--threads Number of threads for parallel processing
--readable_json Save JSON in a human-readable format (True/False, default: False)
"""

from collections import defaultdict
import argparse
from concurrent.futures import ThreadPoolExecutor
from io_handler import MSpangepopDataHandler, MScompute, MSsuccess, MSerror

class TreeNode:
    """Represents a node in the binary tree."""
    
    def __init__(self, node_id, time, interval):
        self.id = node_id  # Node ID
        self.time = time  # Time attribute for ordering nodes
        self.left = None
        self.right = None
        self.mutations = []  # Store mutation details
        self.interval = interval  # Interval for the node
        self.affected_nodes = set()  # Set of all descendant nodes, including itself

class BinaryTree:
    """Handles tree construction and processing."""

    def __init__(self, tree_data):
        self.nodes = {}  # Dictionary to store TreeNode objects by ID
        self.root = None
        self.tree_index = tree_data.get("tree_index", "unknown")  # Ensure default value if missing
        self.interval = tree_data.get("interval", None)  # Retrieve interval if present
        self.lineages = set()  # This to store lineages
        self.build_tree(tree_data)  # Construct tree from input data
    
    def build_tree(self, tree_data):
        """Builds the binary tree from input JSON data with error handling."""
        try:
            if "nodes" not in tree_data or "edges" not in tree_data or "mutations" not in tree_data:
                raise MSerror(f"Malformed tree data for tree index {self.tree_index}. Missing required fields.")

            if self.interval is None:
                raise MSerror(f"Missing interval for tree index {self.tree_index}.")

            for node in tree_data["nodes"]:
                interval = node.get("interval", self.interval)  # Node specific interval or tree interval
                self.nodes[node["id"]] = TreeNode(node["id"], node["time"], interval)

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

            for mutation in tree_data["mutations"]:
                node_id = mutation["node"]
                if node_id not in self.nodes:
                    raise MSerror(f"Mutation assigned to non-existent node {node_id}.")
                self.nodes[node_id].mutations.append({
                    "type": None,
                    "start": None,
                    "length": None
                })

            all_children = {edge["child"] for edge in tree_data["edges"]}
            root_candidates = set(self.nodes.keys()) - all_children
            if root_candidates:
                self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]
            else:
                raise MSerror(f"Tree index {self.tree_index} has no identifiable root node.")
            
            self.compute_affected_nodes()
        
        except KeyError as e:
            raise MSerror(f"Key error during tree building for tree {self.tree_index}: {e}")
        except Exception as e:
            raise MSerror(f"Unexpected error while processing tree {self.tree_index}: {e}")

    def compute_affected_nodes(self):
            """Computes and stores the affected nodes for each node in the tree."""
            
            def dfs(node):
                if not node:
                    return set()
                
                # Check if this is a leaf node (lineage)
                if node.left is None and node.right is None:
                    self.lineages.add(node.id)  # Add to lineages set
                
                affected_nodes = {node.id}
                if node.left:
                    affected_nodes |= dfs(node.left)
                if node.right:
                    affected_nodes |= dfs(node.right)
                node.affected_nodes = affected_nodes
                return affected_nodes

            if not self.root:
                raise MSerror(f"Cannot compute affected nodes: No root node found for tree {self.tree_index}.")
            dfs(self.root)

    def get_all_nodes(self):
        """Returns all nodes, including mutated and non-mutated nodes."""
        all_nodes = []

        def traverse(node):
            if node:
                all_nodes.append({
                    "node": node.id,
                    "mutations": node.mutations,
                    "interval": node.interval,
                    "affected_nodes": list(node.affected_nodes)
                })
                traverse(node.left)
                traverse(node.right)
        
        traverse(self.root)
        return {
            "tree_index": self.tree_index,
            "initial_tree_interval": self.interval,
            "lineages": sorted(list(self.lineages)),  # Add sorted lineages here
            "nodes": all_nodes
        }

def process_tree(tree_data):
    try:
        tree = BinaryTree(tree_data)
        return tree.get_all_nodes()
    except MSerror as e:
        raise MSerror(f"Error while processing a tree: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serialize all nodes into JSON with mutation details, intervals, and affected nodes.")
    parser.add_argument("--json", required=True, help="Path to the input JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome number")
    parser.add_argument("--output_file", required=True, help="Path to the output JSON file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False],
                        help="Save JSON in a human-readable format (True/False, default: False).")
    args = parser.parse_args()

    try:
        json_data = MSpangepopDataHandler.read_json(args.json)
        if not isinstance(json_data, list):
            raise MSerror(f"Invalid JSON format in {args.json}: Expected a list of tree objects.")

        MScompute(f"Starting tree sequence traversal for {args.chromosome} with {args.threads} threads")

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            results = list(filter(None, executor.map(process_tree, json_data)))

        if results:
            MSpangepopDataHandler.save_json(results, args.output_file, args.readable_json)
            MSsuccess(f"Serialized traversal for chromosome {args.chromosome}")
        else:
            raise MSerror(f"No valid trees processed for chromosome {args.chromosome}. Output file not saved.")
    except FileNotFoundError:
        raise MSerror(f"Input file not found: {args.json}")
    except Exception as e:
        raise MSerror(f"Critical error during execution: {e}")
