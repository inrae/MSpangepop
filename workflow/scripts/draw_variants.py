"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Augment JSON file with variant type and size..

--json Path to the JSON file containing tree and mutation data
--output Path to the output JSON file where augmented data will be saved
--yaml Path to the YAML configuration file with variant probabilities.
--chromosome Chromosome name
--threads Number of threads for parallel processing
--minimal_variant_size Minimal size for variants generated
--readable_json Save JSON in a human-readable format (True/False, default: False)
"""

import argparse
import random
import os
import time
import concurrent.futures
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning

import random

class Selector:
    """Static methods for selecting mutation attributes."""
    
    @staticmethod
    def position(interval) -> int:
        """Returns a random position within the given interval."""
        return random.randint(interval[0] , interval[1]-2)

    @staticmethod
    def type(variant_probabilities) -> str:
        """Selects a variant type based on predefined probabilities."""
        try:
            return random.choices(list(variant_probabilities.keys()), weights=variant_probabilities.values())[0]
        except Exception as e:
            raise MSerror(f"Error in selecting type length: {e}")

    @staticmethod 
    def length(length_df, max_length, minimal_variant_size):
        try:
            rand_val = random.random()
            row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]
            lower_bound, upper_bound = map(float, row['size_interval'].strip('[]').split(','))
            length = random.randint(int(lower_bound), int(upper_bound))
            if max_length:
                return max(min(minimal_variant_size, max_length), min(length, max_length))
            else :
                return max(minimal_variant_size, length)
        except Exception as e:
            raise MSerror(f"Error in selecting variant length : {e}")
        
def augment_type(tree, variant_probabilities):
    """Augments mutations in a tree with variant type"""
    try:
        for node in tree.get("nodes", []):
            for mutation in node.get("mutations", []):
                # Assign variant type
                mutation["type"] = Selector.type(variant_probabilities)
        return tree
    except Exception as e:
        MSerror(f"Error augmenting mutations: {e}")
        return None
    
def augment_length_and_position(tree, length_files, minimal_variant_size):
    """Augments mutations in a tree with length and position based on variant types."""
    try:
        for node in tree.get("nodes", []):  
            for mutation in node.get("mutations", []):  # Iterate over mutations of the node

                # Select a random position within the node's interval
                start_pos = Selector.position(node["interval"])

                # Select the mutation length
                variant_type = mutation["type"]

                if variant_type == "SNP": #For SNPs do nothing
                    length = None 

                if variant_type == "INS" : #For INSs select a random length (no limit)
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, None, minimal_variant_size)

                    affected_nodes = node.get("affected_nodes", [])
                    for affected_node_id in affected_nodes:
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] += length #Expend affected node length
                    
                if variant_type == "DUP" : #For INSs select a length with the remaining 3' bases 
                    length_df = length_files.get(variant_type)
                    max_size = node["interval"][1] - start_pos - 1
                    length = Selector.length(length_df, max_size, minimal_variant_size)
                    affected_nodes = node.get("affected_nodes", [])
                    for affected_node_id in affected_nodes:
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] += length
                
                elif variant_type == "DEL":
                    current_size = node["interval"][1] - node["interval"][0]
                    
                    if current_size > 3:  # Ensure the sequence is never less than 2 bases
                        length_df = length_files.get(variant_type)
                        max_size = node["interval"][1] - start_pos - 1
                        length = Selector.length(length_df, max_size, 0)

                        if current_size - length >= 2:  # Ensure deletion does not reduce size below 2
                            affected_nodes = node.get("affected_nodes", [])
                            for affected_node_id in affected_nodes:
                                affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                                if affected_node:
                                    affected_node["interval"][1] -= length  # Subtract for DELs
                        else:
                            length = None  # Skip deletion if it would make sequence < 2
                    else:
                        length = None  # Skip deletion if there are only 2 bases
                        
                if variant_type == "INV" : #For the INVs we dont change the total size
                    length_df = length_files.get(variant_type)
                    max_size = node["interval"][1] - start_pos - 1 
                    length = Selector.length(length_df, max_size, minimal_variant_size)

                mutation["length"] = length
                mutation["start"] = start_pos

        return tree
    except Exception as e:
        raise MSerror(f"Error augmenting mutations with length and position: {e}")

def process_single_tree(tree, length_files, variant_probabilities, minimal_variant_size):
    """Processes a single tree, augmenting mutation types, positions, and lengths."""
    try:
        # First augment with mutation types
        augmented_tree = augment_type(tree, variant_probabilities)
        
        # Now augment with lengths and positions
        augmented_tree = augment_length_and_position(augmented_tree, length_files, minimal_variant_size)
        
        return augmented_tree
    except MSerror as e:
        raise MSerror(f"Error processing the tree: {e}")

def process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_variant_size):
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_single_tree, tree_list, [length_files]*len(tree_list), [variant_probabilities]*len(tree_list), [minimal_variant_size]*len(tree_list)))
    return results 

def main(json_file, output_json_file, yaml_file, chromosome, num_threads, minimal_variant_size, readable_json):
    """Processes a JSON list of trees, augmenting mutations with variant type and size."""
    try:
        start_time = time.time()
        MScompute(f"Generating variants for chromosome {chromosome} using {num_threads} threads")
        
        # Read variant probabilities (only once)
        variant_probabilities = MSpangepopDataHandler.read_yaml(yaml_file)
        if not isinstance(variant_probabilities, dict) or not variant_probabilities:
            raise MSerror('Variant probabilities file is empty or improperly formatted.')

        # Read tree data (only once)
        tree_list = MSpangepopDataHandler.read_json(json_file)
        if not isinstance(tree_list, list) or not tree_list:
            raise MSerror('Tree JSON file is empty or improperly formatted.')

        # Read length distribution files (only once)
        length_files = {}
        for var_type, file_path in {
            'DEL': 'simulation_data/test.tsv',
            'INS': 'simulation_data/test.tsv',
            'INV': 'simulation_data/test.tsv',
            'DUP': 'simulation_data/test.tsv'
        }.items():
            if not os.path.exists(file_path):
                raise MSerror("Length distribution file missing for {var_type}.")
            length_files[var_type] = MSpangepopDataHandler.read_variant_length_file(file_path)

        # Process and augment mutations in parallel
        tree_list = process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_variant_size)

        # Count total mutations
        total_mutations = sum(len(node.get("mutations", [])) for tree in tree_list for node in tree.get("nodes", []))
        
        # Save output
        MSpangepopDataHandler.save_json(tree_list, output_json_file, readable_json=readable_json)

        # Print summary
        end_time = time.time()
        elapsed_time = end_time - start_time
        MSsuccess(f"Successfully processed chromosome {chromosome}, {total_mutations} mutations handled in {elapsed_time:.2f} sec.")

    except Exception as e:
        raise MSerror(f"Critical error processing (Chromosome {chromosome}): {e}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--yaml", required=True, help="Path to the YAML configuration file with variant probabilities.")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--minimal_variant_size", type=int, default=1, help="Minimal size for variants generated")
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False],
                        help="Save JSON in a human-readable format (True/False, default: False).")

    args = parser.parse_args()

    if args.minimal_variant_size < 1:
            raise MSerror('minimal_variant_size cant be less than 1')

    main(args.json, args.output, args.yaml, args.chromosome, args.threads, args.minimal_variant_size, args.readable_json)
