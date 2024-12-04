import argparse
import random
from readfile import read_json, save_json, read_variant_length_file, read_yaml

def set_length(length_df, max_length):
    """
    Determines the variant length based on the length distribution provided and tree boundaries.
    """
    rand_val = random.random()
    row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]

    interval = row['size_interval']
    lower_bound, upper_bound = map(float, interval.strip('[]').split(','))

    length = random.randint(int(lower_bound), int(upper_bound))

    # Ensure the length does not exceed the tree boundary
    return min(length, max_length)


def select_variant_type(variant_probabilities):
    """
    Selects a variant type based on predefined probabilities.
    """
    return random.choices(list(variant_probabilities.keys()), weights=variant_probabilities.values())[0]


def augment_tree_mutations(tree, length_files, variant_probabilities):
    """
    Augments mutations in a single tree with variant type and length based on distributions,
    ensuring the variant length stays within the tree boundaries.
    """
    for mutation in tree["mutations"]:
        # Select the variant type (SNP, Deletion, Insertion)
        variant_type = select_variant_type(variant_probabilities)
        mutation["variant_type"] = variant_type  # Add the variant type to the mutation

        # SNPs have a fixed length of 1
        if variant_type == "SNP":
            mutation["SV_len"] = 1
        else:
            # Retrieve length distribution for the variant type
            length_df = length_files.get(variant_type)
            if length_df is None:
                raise ValueError(f"No length distribution data available for variant type '{variant_type}'.")

            # Determine the maximum allowable length within tree boundaries
            max_length = tree["interval"][1] - mutation["site_position"]

            # Set the length based on the distribution and tree boundaries
            mutation["SV_len"] = set_length(length_df, max_length)
    return tree


def main(json_file, output_json_file, yaml_file):
    """
    Main function that processes a JSON list of trees,
    augmenting mutations in each tree with variant type and size.
    """
    # Read the variant probabilities from the YAML file
    variant_probabilities = read_yaml(yaml_file)

    # Read the JSON file containing tree and mutation data
    tree_list = read_json(json_file)  # This is now a list of trees

    # Read the length distribution files for Deletions and Insertions
    length_files = {
        'Deletion': read_variant_length_file('simulation_data/size_distribDEL.tsv'),
        'Insertion': read_variant_length_file('simulation_data/size_distribINS.tsv'),
    }

    # Process and augment mutations in each tree
    for tree in tree_list:
        augment_tree_mutations(tree, length_files, variant_probabilities)

    # Save the augmented JSON data to the specified output file
    save_json(tree_list, output_json_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment the JSON file with variant type and size.")
    parser.add_argument("--json", help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--yaml", help="Path to the YAML configuration file with variant probabilities.")
    
    args = parser.parse_args()
    main(args.json, args.output, args.yaml)
