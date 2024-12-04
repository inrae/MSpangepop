import argparse
import json
import random
import yaml

from readfile import read_variant_length_file, read_yaml

def read_fai_to_dict(fai_path):
    """
    Reads an FAI file into a dictionary.
    
    Parameters:
        fai_path (str): Path to the FAI file.
        
    Returns:
        dict: A dictionary where the keys are chromosome names (CHROM),
              and the values are dictionaries containing LENGTH, OFFSET,
              LINEBASES, and LINEWIDTH.
    """
    names = ["CHROM", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"]
    fai_dict = {}
    
    with open(fai_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            chrom_info = dict(zip(names[1:], map(int, fields[1:])))
            fai_dict[fields[0]] = chrom_info
            
    return fai_dict

def get_chromosome_boundaries(fai_dict, chromosome):
    """
    Returns the boundaries of a specified chromosome from the FAI dictionary.
    
    Parameters:
        fai_dict (dict): Dictionary containing chromosome information.
        chromosome (str): The chromosome name to look up.
        
    Returns:
        tuple: A tuple containing the start (OFFSET) and end (OFFSET + LENGTH) positions.
               Returns None if the chromosome is not found.
    """
    chrom_data = fai_dict.get(chromosome)
    if chrom_data is None:
        return None  # Chromosome not found
    
    start = chrom_data["OFFSET"]
    end = chrom_data["OFFSET"] + chrom_data["LENGTH"]
    return start, end


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

def assign_mutation_type(yaml_data):
    """
    Selects a mutation type based on the probabilities in the YAML data.
    
    Parameters:
        yaml_data (dict): The data read from the YAML file, containing the variant percentages.
        
    Returns:
        str: The type of mutation.
    """
    # Create a list of mutation types and their corresponding probabilities
    mutation_types = list(yaml_data.keys())
    probabilities = list(yaml_data.values())
    
    # Randomly select a mutation type based on the given probabilities
    mutation_type = random.choices(mutation_types, probabilities)[0]
    
    return mutation_type


def main(json_file, fai_file, output_json_file, yaml_file):
    # Read the FAI file into a dictionary
    fai_data = read_fai_to_dict(fai_file)

    # Read the JSON file
    parsed_data = read_json_file(json_file)

    # Read the YAML file to get mutation percentages
    yaml_data = read_yaml(yaml_file)

    # Extract the chromosome name from the JSON
    chromosome_name = parsed_data["chromosome"]

    # Get the boundaries of the chromosome using the chromosome name
    boundaries = get_chromosome_boundaries(fai_data, chromosome_name)

    # Print the boundaries
    if boundaries:
        print(f"Boundaries for {chromosome_name}: Start = {boundaries[0]}, End = {boundaries[1]}")
    else:
        print(f"Chromosome {chromosome_name} not found in FAI data.")

    # Assign mutation types to each mutation
    for mutation in parsed_data["mutations"]:
        mutation_type = assign_mutation_type(yaml_data)
        mutation["mutation_type"] = mutation_type  # Add the mutation type to each mutation

    # Save the JSON data to the specified output file
    save_json_file(parsed_data, output_json_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment the JSON file with variant type")
    parser.add_argument("--json", help="Path to the JSON file.")
    parser.add_argument("--fai", help="Path to the FAI file.")
    parser.add_argument("--output", help="Path to the output file")
    parser.add_argument("--yaml", help="Path to the YAML file containing variant type percentages")
    args = parser.parse_args()
    main(args.json, args.fai, args.output, args.yaml)
