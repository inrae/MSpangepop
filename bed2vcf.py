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

# génère une séquence pour insertion
def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

def SNP(aa):
	s = 'ACGT'
	s = s.replace(aa, '')
	return(random.choice(s))

# check if reverse
def is_reverse(s):
	return(s=="reverse")

# reverse sequence
def reverse(s):
	seq = Seq(str(s))
	return(seq.reverse_complement())

# create deletion by switching alt and ref seq
def deletion(ref_seq, alt_seq):
	ref = str(ref_seq) + str(alt_seq)
	alt = str(ref_seq)
	return(ref,alt)

def set_ref_alt(ref, alt, row, df):
	df.at[row, "REF"] = ref
	df.at[row, "ALT"] = alt

def get_random_len(svtype):	
	# INS, DEL, DUP, CNV, INV
	df = pd.read_csv("sv_distributions/size_distrib" + svtype + ".tsv", sep="\t")
	pb = df["pb"].tolist()
	li = np.random.multinomial(1, pb)
	i = np.argmax(li)
	interval = df["size_interval"].iloc[i]
	interval = json.loads(interval)
	s = np.random.uniform(interval[0], interval[1], 1).round()
	return(int(s[0]))

### get the sequence of each variants in BED, output VCF
# vcf_df : msprime VCF
# bed_def : VISOR BED
# fa_dict : FASTA where variants will be generated
def get_seq(vcf_df, bed_df, fa_dict, output_file):
	cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
	ncol = len(vcf_df.columns) - len(cols)
	li = [i for i in range(ncol)]
	cols.extend(li)
	vcf_df.columns = cols

	df = vcf_df.join(bed_df)
	# number of variants
	n = len(df)
	# add variants
	series = []

	for i in range(n):
		## nom du chr pris dans le VCF
		chr_name = df["CHROM"].iloc[i]
		## position du variant pris dans le VCF
		start = df["POS"].iloc[i]
		start = int(start)-1 # pour ajuster à l'index python
		## type de SV à générer
		sv_info = df['SVINFO'].iloc[i]
		sv_type = df['SVTYPE'].iloc[i]

		fasta_seq = fa_dict[chr_name].upper()

		if sv_type == "SNP":
			ref = str(fasta_seq.seq[start])
			alt = SNP(ref)

		elif sv_type == "deletion":
			end = start + get_random_len("DEL")
			ref = str(fasta_seq.seq[start:end])
			alt = str(fasta_seq.seq[start])

		elif sv_type == "insertion":
			end = start + get_random_len("INS")
			alt_seq = DNA(end)
			ref = str(fasta_seq.seq[start])
			alt = ref + str(alt_seq)

		elif sv_type == "inversion":
			end = start + get_random_len("INV")
			ref = str(fasta_seq.seq[start:end])
			alt = str(fasta_seq.seq[start]) + reverse(str(fasta_seq.seq[start+1:end]))

		elif re.search("duplication", sv_type):
			end = start + get_random_len("DUP")
			ref = str(fasta_seq.seq[start])
			alt_seq1 = str(fasta_seq.seq[start+1:end])
			cp = sv_info-1
			# multiplication pour obtenir les copies
			alt_seq = alt_seq1*cp
			if sv_type == "inverted tandem duplication":
				alt_seq = reverse(alt_seq)
			
			alt = ref + alt_seq

		elif re.search("translocation", sv_type):
			# TODO choisir la distribution utilisée
			l = get_random_len("INS")
			end = start + l
			
			# informations sur la translocation
			trans_chr = sv_info[0]
			trans_start = sv_info[1]-1
			trans_end = trans_start + l
			fasta_seq_trans = fa_dict[trans_chr].upper()
			
			# translocation réciproque
			if sv_type == "reciprocal translocation":
				ref = str(fasta_seq.seq[start:end])
				ref2 = str(fasta_seq_trans.seq[trans_start:trans_end])
				
				if sv_info[2] == "reverse":
					alt = str(fasta_seq.seq[start]) + reverse(str(fasta_seq_trans.seq[trans_start+1:trans_end]))
				else :
					alt = str(fasta_seq.seq[start]) + str(fasta_seq_trans.seq[trans_start+1:trans_end])

				if sv_info[3] == "reverse":
					alt2 = str(fasta_seq_trans.seq[trans_start]) + reverse(str(fasta_seq.seq[start+1:end]))
				else:
					alt2 = str(fasta_seq_trans.seq[trans_start]) + str(fasta_seq.seq[start+1:end])
			
			# couper-coller
			elif sv_type == "translocation cut-paste":
				# délétion
				ref = str(fasta_seq.seq[start:end])
				alt = str(fasta_seq.seq[start])
				
				# insertion
				ref2 = str(fasta_seq_trans.seq[trans_start])
				if sv_info[2] == "reverse":
					alt2 = ref2 + reverse(str(fasta_seq.seq[start+1:end]))
				else :
					alt2 = ref2 + str(fasta_seq.seq[start+1:end])

			# copier-coller
			else:
				ref = str(fasta_seq.seq[start])
				alt = SNP(ref)

				ref2 = str(fasta_seq_trans.seq[trans_start])
				alt2 = ref2 + str(fasta_seq.seq[start+1:end])
			
			new_vcf_var = vcf_df.iloc[i].tolist()
			new_vcf_var[0] = trans_chr
			new_vcf_var[1] = trans_start+1
			new_vcf_var[3] = ref2
			new_vcf_var[4] = alt2
			save_series = pd.Series(new_vcf_var, index=cols)
			series.append(save_series)

		set_ref_alt(ref, alt, i, df)

	df_add = pd.DataFrame(series)
	df = pd.concat([df, df_add], ignore_index=True)
	df = df.drop(['SVTYPE', 'SVINFO'], axis=1)

	# remove unnecessary columns for VCF
	df["ID"] = "."
	df["QUAL"] = 99
	# output VCF
	if not os.path.exists("results"):
		os.mkdir("results")
	df.to_csv("results/" + output_file + ".vcf", sep="\t", header=False, index=False)