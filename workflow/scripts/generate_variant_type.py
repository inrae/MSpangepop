import argparse
import json
import random
import pandas as pd
import yaml

def read_fai(fai_file):
    """Read a .fai file and return a dictionary of chromosome lengths."""
    try:
        fai_df = pd.read_table(fai_file, header=None, names=["CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"])
        return dict(zip(fai_df["CHROM"], fai_df["LENGTH"]))
    
    except Exception as e:
        print(f"Error reading FAI file: {e}")
        raise

def get_chromosome_lenght(fai_file, chromosome):
    chrom_lengths =read_fai(fai_file)
    return chrom_lengths[chromosome]

def read_json_file(json_path):
    """
    Reads a JSON file containing chromosome, nodes, edges, and mutations information.
    
    Parameters:
        json_path (str): Path to the JSON file.
        
    Returns:
        dict: A dictionary representation of the JSON file.
    """
    with open(json_path, 'r') as file:
        data = json.load(file)
    return data

def save_json_file(data, output_path):
    """
    Saves the data as a JSON file to the specified output path.
    
    Parameters:
        data (dict): The data to save to the file.
        output_path (str): The path where the file should be saved.
    """
    with open(output_path, 'w') as file:
        json.dump(data, file, indent=4)
    print(f"JSON data has been saved to {output_path}")

def read_variant_length_file(file_path):
    """Read length distribution file and parse intervals with probabilities."""
    try:
        df = pd.read_table(file_path)
        df['cumulative_pb'] = df['pb'].cumsum()
        return df
    
    except Exception as e:
        print(f"Error reading variant length file {file_path}: {e}")
        raise

def set_length(length_df, max_length):
    """Determine variant length based on provided length distribution."""

    rand_val = random.random()
    row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]

    interval = row['size_interval']

    lower_bound, upper_bound = map(float, interval.strip('[]').split(','))

    length = random.randint(int(lower_bound), int(upper_bound))

    return min(length, max_length)


def select_variant_type(variant_probabilities):
    """Select a variant type based on predefined probabilities."""
    return random.choices(list(variant_probabilities.keys()), weights=list(variant_probabilities.values()))[0]

def augment_mutations(mutations, length_files, max_length, variant_probabilities):
    """Augment mutations with a randomly selected variant type and length based on probability distributions."""
    for mutation in mutations:
        # Select the variant type (SNP, Deletion, Insertion)
        variant_type = select_variant_type(variant_probabilities)
        mutation["variant_type"] = variant_type  # Add the variant type to the mutation
        
        # Select the length based on the variant type
        if variant_type == "SNP":
            mutation["SV_len"] = 1  # SNP has length 1
        else:
            length_df = length_files.get(variant_type)
            mutation["SV_len"] = set_length(length_df,max_length)
    return mutations

def probabilities_from_yaml(yaml_file):
    """Reads variant probabilities from a YAML configuration file."""
    try:
        with open(yaml_file, 'r') as file:
            variant_probabilities = yaml.safe_load(file)
        
        # Ensure probabilities sum to 100
        if sum(variant_probabilities.values()) != 100:
            raise ValueError("Sum of variant probabilities in YAML must equal 100.")
        
        return variant_probabilities
    except Exception as e:
        print(f"Error reading YAML file: {e}")
        raise

def main(json_file, fai_file, output_json_file, yaml_file):
    # Read the variant probabilities from the YAML file
    variant_probabilities = probabilities_from_yaml(yaml_file)

    # Read the JSON file
    parsed_data = read_json_file(json_file)

    # Extract the chromosome name from the JSON
    chromosome_name = parsed_data["chromosome"]
    max_length = get_chromosome_lenght(fai_file, chromosome_name)

    # Read the length distribution files for Deletions and Insertions
    length_files = {
        'Deletion': read_variant_length_file('simulation_data/size_distribDEL.tsv'),
        'Insertion': read_variant_length_file('simulation_data/size_distribINS.tsv'),
    }

    # Augment mutations with variant type and length
    parsed_data['mutations'] = augment_mutations(parsed_data['mutations'], length_files, max_length, variant_probabilities)

    # Save the JSON data to the specified output file
    save_json_file(parsed_data, output_json_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment the JSON file with variant type and size.")
    parser.add_argument("--json", help="Path to the JSON file.")
    parser.add_argument("--fai", help="Path to the FAI file.")
    parser.add_argument("--output", help="Path to the output JSON file.")
    parser.add_argument("--yaml", help="Path to the YAML configuration file.")
    args = parser.parse_args()
    main(args.json, args.fai, args.output, args.yaml)
