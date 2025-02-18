"""
Author: Lucien Piat
Creation: 4 Feb 2025
Updated: 12 Feb 2025
Institution: INRAe
Project: PangenOak
"""

import argparse
import random
import os
import sys
import time
import concurrent.futures
from readfile import read_json, save_json, read_variant_length_file, read_yaml

def set_length(length_df, max_length, minimal_variant_size):
    """Determines the variant length based on the length distribution provided and tree boundaries."""
    try:
        rand_val = random.random()
        row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]
        lower_bound, upper_bound = map(float, row['size_interval'].strip('[]').split(','))
        length = random.randint(int(lower_bound), int(upper_bound))
        return max(min(minimal_variant_size, max_length), min(length, max_length))  # Ensure length does not exceed the boundary and is not 0
    except Exception as e:
        print(f"❌ MSpangepop -> Error in selecting variant length : {e}", file=sys.stderr)
        sys.exit(1)

def select_variant_type(variant_probabilities):
    """Selects a variant type based on predefined probabilities."""
    try:
        return random.choices(list(variant_probabilities.keys()), weights=variant_probabilities.values())[0]
    except Exception as e:
        print(f"❌ MSpangepop -> Error selecting variant type : {e}", file=sys.stderr)
        sys.exit(1)

def augment_tree_mutations(tree, length_files, variant_probabilities, minimal_variant_size):
    """Augments mutations in a tree with variant type and length based on distributions."""
    try:
        for mutation in tree.get("mutations", []):
            variant_type = select_variant_type(variant_probabilities)
            mutation["variant_type"] = variant_type

            if variant_type == "SNP":
                mutation["len"] = 1
            else:
                length_df = length_files.get(variant_type)
                if length_df is None:
                    raise ValueError(f"❌ MSpangepop -> No length distribution data for variant type '{variant_type}'.")
                max_length = tree["interval"][1] - mutation["site_position"]
                mutation["len"] = set_length(length_df, max_length, minimal_variant_size)

        return tree

    except Exception as e:
        print(f"❌ MSpangepop -> Error augmenting mutations : {e}", file=sys.stderr)
        return None

def process_single_tree(tree, length_files, variant_probabilities, minimal_variant_size):
    return augment_tree_mutations(tree, length_files, variant_probabilities, minimal_variant_size)

def process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_variant_size):
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_single_tree, tree_list, [length_files]*len(tree_list), [variant_probabilities]*len(tree_list), [minimal_variant_size]*len(tree_list)))
    return [tree for tree in results if tree is not None]  # Filter out None results

def main(json_file, output_json_file, yaml_file, chromosome, num_threads, minimal_variant_size):
    """Processes a JSON list of trees, augmenting mutations with variant type and size."""
    try:
        start_time = time.time()
        print(f"🔹 MSpangepop -> Generating variants for chromosome {chromosome} using {num_threads} threads")
        
        if minimal_variant_size < 1:
            raise ValueError(f"\n❌ MSpangepop -> minimal_variant_size cant be less than 1")
        if not os.path.exists(json_file):
            raise FileNotFoundError(f"\n❌ MSpangepop -> Input JSON file not found: {json_file}")
        if not os.path.exists(yaml_file):
            raise FileNotFoundError(f"\n❌ MSpangepop -> YAML configuration file not found: {yaml_file}")

        # Read variant probabilities (only once)
        variant_probabilities = read_yaml(yaml_file)
        if not isinstance(variant_probabilities, dict) or not variant_probabilities:
            raise ValueError("\n❌ MSpangepop -> Variant probabilities file is empty or improperly formatted.")

        # Read tree data (only once)
        tree_list = read_json(json_file)
        if not isinstance(tree_list, list) or not tree_list:
            raise ValueError("\n❌ MSpangepop -> Tree JSON file is empty or improperly formatted.")

        # Read length distribution files (only once)
        length_files = {}
        for var_type, file_path in {
            'DEL': 'simulation_data/test.tsv',
            'INS': 'simulation_data/test.tsv',
            'INV': 'simulation_data/test.tsv',
            'DUP': 'simulation_data/test.tsv'
        }.items():
            if not os.path.exists(file_path):
                print(f"⚠️ MSpangepop -> Length distribution file missing for {var_type}.", file=sys.stderr)
                continue
            length_files[var_type] = read_variant_length_file(file_path)

        # Process and augment mutations in parallel
        tree_list = process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_variant_size)

        # Count total mutations
        total_mutations = sum(len(tree.get("mutations", [])) for tree in tree_list)

        # Save output
        save_json(tree_list, output_json_file)

        # Print summary
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"✅ MSpangepop -> Successfully processed chromosome {chromosome}, {total_mutations} mutations handled in {elapsed_time:.2f} sec.")

    except Exception as e:
        print(f"❌ MSpangepop -> Critical error processing (Chromosome {chromosome}): {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--yaml", required=True, help="Path to the YAML configuration file with variant probabilities.")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--minimal_variant_size", type=int, default=1, help="Minimal size for variants generated")

    args = parser.parse_args()
    main(args.json, args.output, args.yaml, args.chromosome, args.threads, args.minimal_variant_size)
