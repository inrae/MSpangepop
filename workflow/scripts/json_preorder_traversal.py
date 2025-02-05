import argparse
from collections import defaultdict

from readfile import read_json

class TreeNode:
    """Represents a node in the binary tree with an ID, time, children, and mutations."""
    def __init__(self, node_id, time):
        self.id = node_id
        self.time = time
        self.left = None  # Left child (older child)
        self.right = None  # Right child (younger child)
        self.mutations = []  # Mutations on this node

    def __repr__(self):
        return f"TreeNode(id={self.id}, time={self.time}, mutations={self.mutations})"
    
class BinaryTree:
    """Builds and processes a binary tree based on given nodes, edges, and mutations."""
    
    def __init__(self, tree_data):
        self.nodes = {}  # Maps node_id to TreeNode
        self.root = None
        self.build_tree(tree_data)
    
    def build_tree(self, tree_data):
        """Constructs the binary tree from JSON data."""
        # Create nodes
        for node in tree_data["nodes"]:
            self.nodes[node["id"]] = TreeNode(node["id"], node["time"])

        # Build edges (assign left/right children)
        child_counts = defaultdict(int)  # Track children for left/right assignment
        for edge in tree_data["edges"]:
            parent_id, child_id = edge["parent"], edge["child"]
            parent = self.nodes[parent_id]
            child = self.nodes[child_id]

            if child_counts[parent_id] == 0:
                parent.left = child  # Assign first child to left
            else:
                parent.right = child  # Assign second child to right
            child_counts[parent_id] += 1

        # Attach mutations to nodes
        for mutation in tree_data["mutations"]:
            node_id = mutation["node"]
            self.nodes[node_id].mutations.append(mutation)

        # Find the root (oldest node)
        all_children = {edge["child"] for edge in tree_data["edges"]}
        root_candidates = set(self.nodes.keys()) - all_children
        if root_candidates:
            self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]  # Oldest root

    def preorder_traversal(self, node):
        """Performs preorder traversal and collects mutations."""
        if node:
            # Visit the current node first
            print(f"\tNode {node.id} - Mutations: {len(node.mutations)}")
            mutations = node.mutations
            for mutation in mutations: 
                print(f"\t\t - {mutation['variant_type']} at position {mutation['site_position']} at time {mutation['time']}")

            # Traverse the left subtree
            self.preorder_traversal(node.left)
            # Traverse the right subtree
            self.preorder_traversal(node.right)

    def print_preorder(self):
        """Prints the preorder traversal and collected mutations from oldest to newest."""
        if not self.root:
            print("❌ MSpangepop -> No tree found!")
            return
    
        self.preorder_traversal(self.root)


def main(json_file, chromosome):
    """Processes all trees in the input JSON file."""
    json_data = read_json(json_file)
    print(f"\n🔹 MSpangepop -> Starting to process preorder traversal for chr {chromosome}")
    for tree_data in json_data:
        print(f"\tMSpangepop -> Tree {tree_data['tree_index']} {tree_data['interval']}")
        tree = BinaryTree(tree_data)
        tree.print_preorder()
        print(f"\t✅ MSpangepop -> Tree {tree_data['tree_index']} traversed sucessfuly")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate the preorder path from the json file")
    parser.add_argument("--json", required=True, help="Path to the json file.")
    parser.add_argument("--chromosome", required=True, help="chr number")
    args = parser.parse_args()

    main(args.json, args.chromosome)
