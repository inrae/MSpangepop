import pandas as pd
import yaml
import random
import re

def read_yaml(file_path):
    """
    Read a YAML file containing variant information and validate the percentages.
    
    Returns:
        dict: Parsed YAML data.
    """
    with open(file_path, 'r') as stream:
        data = yaml.safe_load(stream)
    
	# Raises ValueError: If the sum of variant percentages is not 100.
    if sum(data.values()) != 100:
        raise ValueError("Sum of values must be 100.")
    
    return data

# Predefined information for each variant type
variant_info_dict = {
    'deletion': None,
    'insertion': None,
    'inversion': None,
    'SNP': None,
    'tandem duplication': 2,
    'inverted tandem duplication': 2,
    'translocation copy-paste': 0,
    'translocation cut-paste': 0,
    'reciprocal translocation': 0
}

def select_chr(fai_file):
    """
    Randomly select a chromosome and position from the FAI file.

    Returns:
        list: Selected chromosome and random position on that chromosome.
    """
    df = pd.read_table(fai_file, header=None, usecols=[0, 1], names=["name", "length"])
    chromosome_names = df["name"].to_list()
    
    selected_chr = random.choice(chromosome_names)
    length = df.loc[df['name'] == selected_chr, 'length'].values[0]
    position = random.randrange(length)
    
    return [selected_chr, position]

def select_orientation(info_list):
    """
    Randomly select the orientation (forward or reverse) for a variant.

    Returns:
        list: Updated list including the orientation.
    """
    info_list.append(random.choice(["forward", "reverse"]))
    return info_list

def generate_type(num_variants, yaml_file, fai_file):
    """
    Generate a DataFrame containing information about structural variants.

    Returns:
        pd.DataFrame: DataFrame containing variant types and their information.
    """
    variant_distribution = read_yaml(yaml_file)
    
    variant_types = []
    variant_infos = []
    
    # Filter the dictionary to retain only non-zero values
    filtered_variants = {k: v for k, v in variant_distribution.items() if v != 0}
    
    for variant_type, percentage in filtered_variants.items():
        # Calculate the number of variants to create based on percentage
        num_to_generate = round(percentage * num_variants / 100)
        
        # Generate information for each variant
        if re.search("translocation", variant_type):
            for _ in range(num_to_generate):
                info = select_chr(fai_file)
                info = select_orientation(info)
                if variant_type == "reciprocal translocation":
                    info = select_orientation(info)
                variant_infos.append(info)
                variant_types.append(variant_type)
        else:
            info = variant_info_dict[variant_type]
            variant_infos.extend([info] * num_to_generate)
            variant_types.extend([variant_type] * num_to_generate)
    
    # If the total generated variants are less than required, fill with SNPs
    if len(variant_types) < num_variants:
        additional_variants = num_variants - len(variant_types)
        variant_types.extend(["SNP"] * additional_variants)
        variant_infos.extend([None] * additional_variants)

    # Create a DataFrame from the generated variant types and their info
    variant_df = pd.DataFrame({'SVTYPE': variant_types, 'SVINFO': variant_infos})

    # Shuffle the DataFrame
    variant_df = variant_df.sample(frac=1, ignore_index=True)

    # Save the DataFrame to a TSV file
    variant_df.to_csv("random_var.tsv", sep="\t", index=False)
    
    return variant_df
