import argparse
from Bio import SeqIO

def read_fasta(input_fasta):
    """
    Reads a FASTA file and returns a dictionary of sequences.
    
    Parameters:
        input_fasta (str): Path to the FASTA file.
        
    Returns:
        dict: A dictionary with sequence IDs as keys and sequence data as values.
    """
    try:
        return SeqIO.index(input_fasta, "fasta")
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        raise

def main(fasta_file):
    fast_sequence = read_fasta(fasta_file)
    print(fast_sequence[:10])   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    
    args = parser.parse_args()
    main(args.fasta)
