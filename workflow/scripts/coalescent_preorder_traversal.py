import argparse
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from io_handler import MSpangepopDataHandler, MScompute, MSsuccess

class TreeNode:
    """Represents a node in the binary tree."""
    
    def __init__(self, node_id, time):
        self.id = node_id  # Node ID
        self.time = time  # Time attribute for ordering nodes
        self.left = None
        self.right = None
        self.mutations = []  # List of mutations at this node
        self.leaves = set()  # Set of all descendant leaf nodes, including itself if it is a leaf

    def __repr__(self):
        return f"TreeNode(id={self.id}, time={self.time}, mutations={self.mutations}, leaves={self.leaves})"

class BinaryTree:
    """Handles tree construction and processing."""

    def __init__(self, tree_data):
        self.nodes = {}  # Dictionary to store TreeNode objects by ID
        self.root = None
        self.tree_index = tree_data["tree_index"]  # Index of the tree (useful if multiple trees exist)
        self.build_tree(tree_data)  # Construct tree from input data

    def build_tree(self, tree_data):
        """Builds the binary tree from input JSON data."""

        # Step 1: Create all nodes without linking them
        for node in tree_data["nodes"]:
            self.nodes[node["id"]] = TreeNode(node["id"], node["time"])

        # Step 2: Link nodes using parent-child relationships
        child_counts = defaultdict(int)  # Track how many children each parent has
        for edge in tree_data["edges"]:
            parent_id, child_id = edge["parent"], edge["child"]
            parent = self.nodes[parent_id]
            child = self.nodes[child_id]

            # Assign child nodes to left or right, depending on order
            if child_counts[parent_id] == 0:
                parent.left = child
            else:
                parent.right = child
            child_counts[parent_id] += 1

        # Step 3: Assign mutations to nodes
        for mutation in tree_data["mutations"]:
            node_id = mutation["node"]
            self.nodes[node_id].mutations.append(mutation)

        # Step 4: Identify the root node
        all_children = {edge["child"] for edge in tree_data["edges"]}
        root_candidates = set(self.nodes.keys()) - all_children
        if root_candidates:
            self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]

        # Step 5: Compute affected leaf sets for each node
        self.compute_leaf_sets()

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

        # Start computing leaf sets from the root
        dfs(self.root)

    def get_all_leaves(self):
        """Returns a sorted list of all unique leaf node IDs in the tree."""
        all_leaves = set()
        for node in self.nodes.values():
            if not node.left and not node.right:  # Check if node is a leaf
                all_leaves.add(node.id)
        return sorted(all_leaves)  # Sorting for consistency

    def preorder_traversal(self, node, traversal):
        """Performs a preorder traversal of the tree, collecting mutations and affected leaves."""
        
        if node:
            traversal.append({
                "node": node.id,
                "mutations": [
                    {
                        "site_position": m["site_position"],
                        "length": m["len"],
                        "variant_type": m["variant_type"]
                    }
                    for m in node.mutations
                ],
                "affected_leaves": list(node.leaves)
            })
            # Recursively visit left and right children
            self.preorder_traversal(node.left, traversal)
            self.preorder_traversal(node.right, traversal)

    def get_preorder_traversal(self):
        """Returns the preorder traversal including affected leaves for each node, along with tree leaves."""
        traversal = []
        self.preorder_traversal(self.root, traversal)
        return {
            "tree_index": self.tree_index,
            "traversal": traversal,
            "leaves_present": self.get_all_leaves()  # Include the list of all leaves in this tree
        }

def process_tree(tree_data):
    """
    Processes a single tree and returns the preorder traversal result.
    Includes mutations, affected leaves, and all leaves present in the tree.
    """
    tree = BinaryTree(tree_data)
    return tree.get_preorder_traversal()

def main(json_file, chromosome, output_file, threads, readble_json):
    """
    Main function to process tree data in parallel and save the results to a JSON file.
    
    Parameters:
        json_file (str): Path to the input JSON file containing tree data.
        chromosome (str): Chromosome identifier.
        output_file (str): Path to save the serialized traversal output.
        threads (int): Number of parallel threads to use.
    """
    
    # Load input data
    json_data = MSpangepopDataHandler.read_json(json_file)
    
    MScompute(f"Processing traversal of chromosome {chromosome} with {threads} threads")

    # Process trees in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_tree, json_data))

    # Save output JSON, which now includes all leaves inside each tree
    MSpangepopDataHandler.save_json(list(results), output_file, readble_json)

    MSsuccess(f"Serialized traversal for chromosome {chromosome}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serialize the preorder traversal including affected leaves into JSON.")
    parser.add_argument("--json", required=True, help="Path to the input JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome number")
    parser.add_argument("--output_file", required=True, help="Path to the output JSON file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    parser.add_argument("--readble_json", default=False)
    args = parser.parse_args()

    # Execute the main function
    main(args.json, args.chromosome, args.output_file, args.threads, args.readble_json)
