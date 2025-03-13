import argparse
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from readfile import read_json, save_json

class TreeNode:
    def __init__(self, node_id, time):
        self.id = node_id
        self.time = time
        self.left = None
        self.right = None
        self.mutations = []

    def __repr__(self):
        return f"TreeNode(id={self.id}, time={self.time}, mutations={self.mutations})"
    
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
            if child_counts[parent_id] == 0:
                parent.left = child
            else:
                parent.right = child
            child_counts[parent_id] += 1

        for mutation in tree_data["mutations"]:
            node_id = mutation["node"]
            self.nodes[node_id].mutations.append(mutation)

        all_children = {edge["child"] for edge in tree_data["edges"]}
        root_candidates = set(self.nodes.keys()) - all_children
        if root_candidates:
            self.root = self.nodes[min(root_candidates, key=lambda n: self.nodes[n].time)]

    def preorder_traversal(self, node, traversal):
        if node:
            traversal_entry = {
                "node": node.id,
                "mutations": []
            }
            for mutation in node.mutations:
                traversal_entry["mutations"].append({
                    "site_position": mutation["site_position"],
                    "length": mutation.get("length", None),
                    "variant_type": mutation["variant_type"]
                })
            
            traversal.append(traversal_entry)
            self.preorder_traversal(node.left, traversal)
            self.preorder_traversal(node.right, traversal)

    def get_preorder_traversal(self):
        traversal = []
        self.preorder_traversal(self.root, traversal)
        return {"tree_index": self.tree_index, "traversal": traversal}

def process_tree(tree_data):
    tree = BinaryTree(tree_data)
    return tree.get_preorder_traversal()

def main(json_file, chromosome, output_file, threads):
    json_data = read_json(json_file)
    
    print(f"🔹Asm4pg -> Processing traversal of chromosome {chromosome} with {threads} threads")
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(process_tree, json_data))
    
    save_json(results, output_file)

    print(f"✅ Asm4pg -> Serialized traversal of chromosome {chromosome}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serialize the preorder traversal of the tree into a JSON file")
    parser.add_argument("--json", required=True, help="Path to the input JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome number")
    parser.add_argument("--output", required=True, help="Path to the output JSON file.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    args = parser.parse_args()

    main(args.json, args.chromosome, args.output, args.threads)
