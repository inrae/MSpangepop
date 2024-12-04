import argparse
import random
from readfile import read_fai, read_json, save_json, read_variant_length_file, read_yaml

def get_chromosome_length(fai_file, chromosome):
    """
    Retrieves the length of a specific chromosome from the FAI file.
    
    Parameters:
        fai_file (str): Path to the FAI file.
        chromosome (str): Chromosome name to get the length for.
    
    Returns:
        int: The length of the specified chromosome.
    
    Raises:
        ValueError: If the specified chromosome is not found in the FAI file.
    """
    chrom_lengths = read_fai(fai_file)
    if chromosome not in chrom_lengths:
        raise ValueError(f"Chromosome '{chromosome}' not found in the FAI file.")
    
    return chrom_lengths[chromosome]


def set_length(length_df, max_length):
    """
    Determines the variant length based on the length distribution provided.
    
    Parameters:
        length_df (pd.DataFrame): DataFrame with variant lengths and probabilities.
        max_length (int): The maximum allowable variant length (chromosome length).
    
    Returns:
        int: A valid variant length that does not exceed max_length.
    """
    rand_val = random.random()
    row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]

    interval = row['size_interval']
    lower_bound, upper_bound = map(float, interval.strip('[]').split(','))

    length = random.randint(int(lower_bound), int(upper_bound))

    return min(length, max_length)

def select_variant_type(variant_probabilities):
    """
    Selects a variant type based on predefined probabilities.
    
    Parameters:
        variant_probabilities (dict): A dictionary with variant types as keys and their probabilities as values.
    
    Returns:
        str: The selected variant type (SNP, Deletion, Insertion, etc.).
    """
    return random.choices(list(variant_probabilities.keys()), weights=variant_probabilities.values())[0]

def augment_mutations(mutations, length_files, max_length, variant_probabilities):
    """
    Augments mutations with a randomly selected variant type and length based on probability distributions.
    
    Parameters:
        mutations (list): List of mutation dictionaries to augment.
        length_files (dict): Dictionary of length distribution DataFrames for variants (Deletion, Insertion).
        max_length (int): The maximum allowable variant length (chromosome length).
        variant_probabilities (dict): Dictionary of variant types and their selection probabilities.
    
    Returns:
        list: The mutated mutations with added variant type and length.
    
    Raises:
        ValueError: If information about a selected variant type is missing.
    """
    for mutation in mutations:
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
            
            # Set the length based on the distribution
            mutation["SV_len"] = set_length(length_df, max_length)

    return mutations

def main(json_file, fai_file, output_json_file, yaml_file):
    """
    Main function that orchestrates the reading of input files, mutation augmentation, and output saving.
    
    Parameters:
        json_file (str): Path to the JSON file.
        fai_file (str): Path to the FAI file.
        output_json_file (str): Path where the augmented JSON data will be saved.
        yaml_file (str): Path to the YAML file with variant probabilities.
    """
    # Read the variant probabilities from the YAML file
    variant_probabilities = read_yaml(yaml_file)

    # Read the JSON file containing mutation data
    parsed_data = read_json(json_file)

    # Extract the chromosome name and its maximum length
    chromosome_name = parsed_data["chromosome"]
    max_length = get_chromosome_length(fai_file, chromosome_name)

    # Read the length distribution files for Deletions and Insertions
    length_files = {
        'Deletion': read_variant_length_file('simulation_data/size_distribDEL.tsv'),
        'Insertion': read_variant_length_file('simulation_data/size_distribINS.tsv'),
    }

    # Augment mutations with variant type and length
    parsed_data['mutations'] = augment_mutations(parsed_data['mutations'], length_files, max_length, variant_probabilities)

    # Save the augmented JSON data to the specified output file
    save_json(parsed_data, output_json_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment the JSON file with variant type and size.")
    parser.add_argument("--json", help="Path to the JSON file containing mutation data.")
    parser.add_argument("--fai", help="Path to the FAI file containing chromosome information.")
    parser.add_argument("--output", help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--yaml", help="Path to the YAML configuration file with variant probabilities.")
    
    args = parser.parse_args()
    main(args.json, args.fai, args.output, args.yaml)
