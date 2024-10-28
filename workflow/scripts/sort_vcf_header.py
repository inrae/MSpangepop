import argparse

def parse_vcf_header(vcf_header):
    """
    Parse the VCF header to extract the general info, contigs, and sort the contigs.

    :param vcf_header: str, the VCF header as a string.
    :return: list, containing combined lines of parsed results.
    """
    # Split the input into lines
    lines = vcf_header.strip().split('\n')

    # Initialize sections
    first_part = []
    contigs = []
    second_part = []
    
    # Iterate through each line and categorize it
    start = True
    for line in lines:
        if line.startswith("##contig"):
            contigs.append(line)
            start = False
        elif start:
            first_part.append(line)
        else:
            second_part.append(line)
    
    def sort_key(contig_line):
        contig_id = contig_line.split('<ID=')[1].split(',')[0]
        if contig_id.startswith("."):
            return (0, contig_id)  # Highest priority
        elif contig_id[0].isdigit():
            return (2, int(contig_id.lstrip('0')))  # Numeric IDs, removing leading zeros for sorting
        else:
            return (1, contig_id)  # Alphabetic IDs

    sorted_contigs = sorted(contigs, key=sort_key)
    
    # Combine all parts into a single list for output
    combined_lines = []
    combined_lines.extend(first_part)
    combined_lines.extend(sorted_contigs)
    combined_lines.extend(second_part)

    return combined_lines

def main(input_file, output_file):
    # Read the VCF header from the input file
    with open(input_file, 'r') as infile:
        vcf_header = infile.read()
        
    # Parse the VCF header
    parsed_result = parse_vcf_header(vcf_header)

    # Write the parsed result line by line to the output file
    with open(output_file, 'w') as outfile:
        for line in parsed_result:
            outfile.write(line + '\n')

if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Parse VCF header and sort contigs.")
    parser.add_argument("input_file", help="Path to the input VCF header file.")
    parser.add_argument("output_file", help="Path to the output file for results.")

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the main function with the input and output file paths
    main(args.input_file, args.output_file)
