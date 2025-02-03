import argparse
import random
import os
import sys
import time
from readfile import read_json, save_json, read_variant_length_file, read_yaml

def set_length(length_df, max_length):
    """Determines the variant length based on the length distribution provided and tree boundaries."""
    try:
        rand_val = random.random()
        row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]

        interval = row['size_interval']
        lower_bound, upper_bound = map(float, interval.strip('[]').split(','))

        length = random.randint(int(lower_bound), int(upper_bound))

        # Ensure the length does not exceed the tree boundary
        return min(length, max_length)
    
    except Exception as e:
        print(f"MSpangepop -> Error in selecting variant length : {e}", file=sys.stderr)
        sys.exit(1)

def select_variant_type(variant_probabilities):
    """Selects a variant type based on predefined probabilities."""
    try:
        return random.choices(list(variant_probabilities.keys()), weights=variant_probabilities.values())[0]
    except Exception as e:
        print(f"MSpangepop -> Error selecting variant type : {e}", file=sys.stderr)
        sys.exit(1)

def augment_tree_mutations(tree, length_files, variant_probabilities):
    """Augments mutations in a tree with variant type and length based on distributions."""
    try:
        for mutation in tree.get("mutations", []):
            # Select the variant type
            variant_type = select_variant_type(variant_probabilities)
            mutation["variant_type"] = variant_type  # Add the variant type to the mutation

            # SNPs have a fixed length of 1
            if variant_type == "SNP":
                mutation["SV_len"] = 1
            else:
                # Retrieve length distribution for the variant type
                length_df = length_files.get(variant_type)
                if length_df is None:
                    raise ValueError(f"MSpangepop -> No length distribution data for variant type '{variant_type}'.")

                # Determine the maximum allowable length within tree boundaries
                max_length = tree["interval"][1] - mutation["site_position"]

                # Set the length based on the distribution and tree boundaries
                mutation["SV_len"] = set_length(length_df, max_length)

        return tree

    except Exception as e:
        print(f"MSpangepop -> Error augmenting mutations : {e}", file=sys.stderr)
        sys.exit(1)

def main(json_file, output_json_file, yaml_file, chromosome):
    """Processes a JSON list of trees, augmenting mutations with variant type and size."""
    try:
        start_time = time.time()
        print(f"MSpangepop -> Generating variants for chromosome {chromosome}")
        if not os.path.exists(json_file):
            raise FileNotFoundError(f"Input JSON file not found: {json_file}")
        if not os.path.exists(yaml_file):
            raise FileNotFoundError(f"YAML configuration file not found: {yaml_file}")

        # Read variant probabilities
        variant_probabilities = read_yaml(yaml_file)
        if not isinstance(variant_probabilities, dict) or not variant_probabilities:
            raise ValueError("Variant probabilities file is empty or improperly formatted.")

        # Read tree data
        tree_list = read_json(json_file)
        if not isinstance(tree_list, list) or not tree_list:
            raise ValueError("Tree JSON file is empty or improperly formatted.")

        # Read length distribution files
        length_files = {}
        for var_type, file_path in {
            'Deletion': 'simulation_data/size_distribDEL.tsv',
            'Insertion': 'simulation_data/size_distribINS.tsv'
        }.items():
            if not os.path.exists(file_path):
                print(f"MSpangepop -> Warning: Length distribution file missing for {var_type}.", file=sys.stderr)
                continue
            length_files[var_type] = read_variant_length_file(file_path)

        # Process and augment mutations in each tree
        for tree in tree_list:
            augment_tree_mutations(tree, length_files, variant_probabilities)

        # Save output
        save_json(tree_list, output_json_file)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"MSpangepop -> Successfully processed chromosome {chromosome} in {elapsed_time/60:.2f} min.. Output saved to {output_json_file}")

    except Exception as e:
        print(f"MSpangepop -> Critical error processing (Chromosome {chromosome}): {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--yaml", required=True, help="Path to the YAML configuration file with variant probabilities.")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    
    args = parser.parse_args()
    main(args.json, args.output, args.yaml, args.chromosome)
