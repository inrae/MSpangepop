import argparse
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from data_handler import MSpangepopDataHandler

class TreeNode:
    def __init__(self, node_id, time):
        self.id = node_id
        self.time = time
        self.left = None
        self.right = None
        self.mutations = []
        self.parents = []  # Track ancestor nodes

    def __repr__(self):
        return f"TreeNode(id={self.id}, time={self.time}, mutations={self.mutations}, parents={self.parents})"
    
class BinaryTree:
    def __init__(self, tree_data):
        self.nodes = {}
        self.root = None
        self.tree_index = tree_data["tree_index"]
        self.build_tree(tree_data)
    
    def build_tree(self, tree_data):
        for node in tree_data["nodes"]:
            self.nodes[node["id"]] = TreeNode(node["id"], node["time"])

        child_counts = defaultdict(int)
        for edge in tree_data["edges"]:
            parent_id, child_id = edge["parent"], edge["child"]
            parent = self.nodes[parent_id]
            child = self.nodes[child_id]

            # Assign left/right children
            if child_counts[parent_id] == 0:
                parent.left = child
            else:
                parent.right = child
            child_counts[parent_id] += 1

            # Track ancestors
            child.parents.append(parent_id)

        for mutation in tree_data["mutations"]:
            node_id = mutation["node"]
            self.nodes[node_id].mutations.append(mutation)

        all_children = {edge["child"] for edge in tree_data["edges"]}
        root_candidates = set(self.nodes.keys()) - all_children
        if root_candidates:
            self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]

    def preorder_traversal(self, node, traversal):
        """Collects preorder traversal with mutation details."""
        if node:
            traversal_entry = {
                "node": node.id,
                "mutations": [{"site_position": m["site_position"], "length": m.get("length"), "variant_type": m["variant_type"]} for m in node.mutations]
            }
            traversal.append(traversal_entry)
            self.preorder_traversal(node.left, traversal)
            self.preorder_traversal(node.right, traversal)

    def get_preorder_traversal(self):
        """Returns preorder traversal of the tree."""
        traversal = []
        self.preorder_traversal(self.root, traversal)
        return {"tree_index": self.tree_index, "traversal": traversal}

    def get_ancestral_data(self):
        """Returns ancestral relationships for all nodes in the tree."""
        return {
            "tree_index": self.tree_index,
            "nodes": [{"node": node_id, "ancestors": node.parents} for node_id, node in self.nodes.items()]
        }

def process_tree(tree_data):
    """Processes a single tree and returns both traversal and ancestral data."""
    tree = BinaryTree(tree_data)
    return tree.get_preorder_traversal(), tree.get_ancestral_data()

def main(json_file, chromosome, traversal_output, ancestry_output, threads):
    json_data = MSpangepopDataHandler.read_json(json_file)
    
    print(f"🔹Asm4pg -> Processing traversal of chromosome {chromosome} with {threads} threads")
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_tree, json_data))
    
    # Separate results into traversal and ancestry data
    traversal_results, ancestry_results = zip(*results)

    MSpangepopDataHandler.save_json(list(traversal_results), traversal_output)
    MSpangepopDataHandler.save_json(list(ancestry_results), ancestry_output)

    print(f"✅ Asm4pg -> Serialized traversal and ancestry of chromosome {chromosome}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serialize the preorder traversal and ancestry information of the tree into JSON files")
    parser.add_argument("--json", required=True, help="Path to the input JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome number")
    parser.add_argument("--traversal_output", required=True, help="Path to the output JSON file for traversal.")
    parser.add_argument("--ancestry_output", required=True, help="Path to the output JSON file for ancestry information.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    args = parser.parse_args()

    main(args.json, args.chromosome, args.traversal_output, args.ancestry_output, args.threads)
