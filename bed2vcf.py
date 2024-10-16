import pandas as pd
import re, random, json, os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

# read FASTA as a BioSeq dictionary
# (not directly loaded in memory)
def read_fa(input_file):
    record_dict = SeqIO.index(input_file, "fasta")
    return record_dict

# Generate a random DNA sequence for insertion
def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

# Generate a SNP
def SNP(aa):
    s = 'ACGT'
    s = s.replace(aa, '')
    return random.choice(s)

# Check if the sequence is reverse
def is_reverse(s):
    return s == "reverse"

# Reverse complement a sequence
def reverse(s):
    seq = Seq(str(s))
    return seq.reverse_complement()

# Create deletion by switching alt and ref sequences
def deletion(ref_seq, alt_seq):
    ref = str(ref_seq) + str(alt_seq)
    alt = str(ref_seq)
    return ref, alt

# Set reference and alternate alleles in the DataFrame
def set_ref_alt(ref, alt, row, df):
    df.at[row, "REF"] = ref
    df.at[row, "ALT"] = alt

# Get random length for structural variants
def get_random_len(svtype):    
    df = pd.read_csv(f"sv_distributions/size_distrib{svtype}.tsv", sep="\t")
    pb = df["pb"].tolist()
    li = np.random.multinomial(1, pb)
    i = np.argmax(li)
    interval = df["size_interval"].iloc[i]
    interval = json.loads(interval)
    s = np.random.uniform(interval[0], interval[1], 1).round()
    return int(s[0])

# Main function to process VCF and BED data, and generate sequence variants
def get_seq(vcf_df, bed_df, fa_dict, output_file):
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    ncol = len(vcf_df.columns) - len(cols)
    li = [i for i in range(ncol)]
    cols.extend(li)
    vcf_df.columns = cols

    df = vcf_df.join(bed_df)
    series = []

    for i in range(len(df)):
        # Get chromosome name, variant position, SV type and sequence
        chr_name = df["CHROM"].iloc[i]
        start = int(df["POS"].iloc[i]) - 1  # Adjust to Python indexing (0-based)
        sv_info = df['SVINFO'].iloc[i]
        sv_type = df['SVTYPE'].iloc[i]

        # Retrieve the reference sequence
        fasta_seq = fa_dict[chr_name].upper()
        fasta_seq_len = len(fasta_seq.seq)

        # Ensure that start and end indices are within bounds
        if start >= fasta_seq_len:
            print(f"Start index {start} is out of bounds for chromosome {chr_name} (length {fasta_seq_len}). Skipping variant.")
            continue

        # Handle different structural variant types
        if sv_type == "SNP":
            ref = str(fasta_seq.seq[start])
            alt = SNP(ref)

        elif sv_type == "deletion":
            end = min(start + get_random_len("DEL"), fasta_seq_len)
            ref = str(fasta_seq.seq[start:end])
            alt = str(fasta_seq.seq[start])

        elif sv_type == "insertion":
            alt_seq = DNA(get_random_len("INS"))
            ref = str(fasta_seq.seq[start])
            alt = ref + alt_seq

        elif sv_type == "inversion":
            end = min(start + get_random_len("INV"), fasta_seq_len)
            ref = str(fasta_seq.seq[start:end])
            alt = str(fasta_seq.seq[start]) + reverse(str(fasta_seq.seq[start+1:end]))

        elif re.search("duplication", sv_type):
            end = min(start + get_random_len("DUP"), fasta_seq_len)
            ref = str(fasta_seq.seq[start])
            alt_seq1 = str(fasta_seq.seq[start+1:end])
            cp = sv_info - 1  # Number of copies
            alt_seq = alt_seq1 * cp
            if sv_type == "inverted tandem duplication":
                alt_seq = reverse(alt_seq)
            alt = ref + alt_seq

        elif re.search("translocation", sv_type):
            l = get_random_len("INS")
            end = min(start + l, fasta_seq_len)

            trans_chr, trans_start, direction1, direction2 = sv_info
            trans_start -= 1
            trans_end = trans_start + l
            fasta_seq_trans = fa_dict[trans_chr].upper()
            trans_seq_len = len(fasta_seq_trans.seq)

            if trans_start >= trans_seq_len or trans_end > trans_seq_len:
                print(f"Translocation start or end out of bounds for chromosome {trans_chr}. Skipping variant.")
                continue

            if sv_type == "reciprocal translocation":
                ref = str(fasta_seq.seq[start:end])
                alt = str(fasta_seq.seq[start]) + reverse(fasta_seq_trans.seq[trans_start+1:trans_end]) if direction1 == "reverse" else str(fasta_seq.seq[start]) + str(fasta_seq_trans.seq[trans_start+1:trans_end])
                ref2 = str(fasta_seq_trans.seq[trans_start])
                alt2 = str(fasta_seq_trans.seq[trans_start]) + reverse(fasta_seq.seq[start+1:end]) if direction2 == "reverse" else str(fasta_seq_trans.seq[trans_start]) + str(fasta_seq.seq[start+1:end])

                # Save the second translocation variant
                new_vcf_var = vcf_df.iloc[i].tolist()
                new_vcf_var[0] = trans_chr
                new_vcf_var[1] = trans_start + 1
                new_vcf_var[3] = ref2
                new_vcf_var[4] = alt2
                series.append(pd.Series(new_vcf_var, index=cols))

            set_ref_alt(ref, alt, i, df)
            continue

        set_ref_alt(ref, alt, i, df)

    df_add = pd.DataFrame(series)
    df = pd.concat([df, df_add], ignore_index=True)
    df = df.drop(['SVTYPE', 'SVINFO'], axis=1)

    # Remove unnecessary columns for VCF and output
    df["ID"] = "."
    df["QUAL"] = 99
    df.to_csv(output_file, sep="\t", header=False, index=False)