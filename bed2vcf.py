import pandas as pd
import re, random, json
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

def set_ref_alt(ref, alt, row, vcf_df):
	vcf_df.at[row, "REF"] = ref
	vcf_df.at[row, "ALT"] = alt

def get_random_len(svtype):
	# INV
	if svtype == "INV":
		rng = np.random.default_rng()
		r = rng.random()
		if r > 0.17:
			s = np.random.uniform(250, 6000,1).round()
		else:
			s = np.random.uniform(6000, 45000,1).round()
	# INS, DEL, DUP, CNV
	else:
		df = pd.read_csv("sv_distributions/size_distrib" + svtype + ".tsv", sep="\t")
		pb = df["pb"].tolist()
		li = np.random.multinomial(1, pb)
		for i in range(len(li)):	
			if li[i] != 0:
				interval = df["size_interval"].iloc[i]
				interval = json.loads(interval)
				s = np.random.uniform(interval[0], interval[1], 1).round()
	return(s[0])

# get the sequence of each variants in BED
# output file is a VCF
# (currently, not all variant supported by VISOR are included)
# TODO : missing VISOR variant type : SNP, MNP and 4 types of tandem repeat

# vcf_df : msprime VCF
# bed_def : VISOR BED
# fa_dict : FASTA where variants will be generated
def get_seq(vcf_df, bed_df, fa_dict, output_file):
	# number of variants
	n = len(vcf_df)

	for i in range(n):
		# pour chaque SV du BED
		sv = bed_df.iloc[i]
		# pour récupérer la séquence
		chr = sv[0]
		start = sv[1]
		end = sv[2]
		# le type de SV
		t = sv[3]
		
		fa_seq = fa_dict[chr]

		alt_seq = fa_seq.seq[start-1:end]
		ref_seq = fa_seq.seq[start-2]
		
		# DEL
		if t == "deletion":
			end = get_random_len("DEL")

			ref = str(ref_seq) + str(alt_seq)
			alt = str(ref_seq)

		# INS
		elif t == "insertion":
			# générer la séquence insérée
			alt_seq = DNA(end-start)
			ref = str(ref_seq)
			alt = str(ref_seq) + str(alt_seq)
			
		# INV
		elif t == "inversion":
			ref = str(ref_seq) + str(alt_seq)
			# reverse alternative sequence
			alt_seq = reverse(alt_seq)
			alt = str(ref_seq) + str(alt_seq)

		# DUP
		# there are 2 duplication types
		elif re.search("duplication", t):
			# get number of copies (minus 1 since it is already in the genome)
			cp = int(sv[4])-1
			# make copy of sequence
			sv_seq = str(alt_seq)*cp
			if t == "inverted tandem duplication":
				# reverse alternative sequence
				sv_seq = reverse(sv_seq)
			ref = str(ref_seq)
			alt = str(ref_seq) + str(sv_seq)
			
		# TRA
		# there are 3 translocation types
		elif re.search("translocation", t):
			# get information field
			infos = sv[4]
			# retrieve each info
			info = infos.split(":")
			tr_start = int(info[2])
			tr_chr = info[1]

			fa_seq2 = fa_dict[tr_chr]

			if info[3] == "reverse":
				# reverse alternative sequence
				alt_seq = reverse(alt_seq)
			
			if t == "translocation cut-paste":
				# cut-paste = INS + DEL
				# DEL at the "cut" (reference) position
				ref = str(ref_seq) + str(alt_seq)
				alt = str(ref_seq)

				# dupliquer la ligne pour ajouter une insertion
				# new row
				new_row = vcf_df.iloc[i].tolist()
				
				# récupère le bon chr
				# translocation start
				
				
				new_row[0] = tr_chr
				
				# translocation position
				new_row[1] = tr_start-2 #???
				# INS at the "paste" position
				ref_seq2 = fa_seq2.seq[tr_start-2]
				alt_seq2 = str(ref_seq2) + str(alt_seq)

				new_row[3] = ref_seq2
				new_row[4] = alt_seq2

				# vcf_df = vcf_df.append(new_row, ignore_index = True)
				vcf_df.loc[len(vcf_df)] = new_row
				# rows_to_add.append(new_sv)
				
			elif t == "reciprocal translocation":
				# à partir de la position de réf
				# translocation 1
				tr_seq = fa_seq.seq[tr_start-1:tr_start+len(alt_seq)]

				ref = str(ref_seq) + str(alt_seq)
				alt = str(ref_seq) + str(tr_seq)
				
				
				# translocation 2
				# à partir de la position de alt
				ref_tr_seq = fa_seq2.seq[tr_start-2]
				ref_seq2 = str(ref_tr_seq) + str(tr_seq)

				if info[4] == "reverse":
					# reverse alternative sequence
					alt_seq = reverse(alt_seq)
	
				alt_seq2 = str(ref_tr_seq) + str(alt_seq)
				# new row
				new_row = vcf_df.iloc[i].tolist()
				new_row[0] = tr_chr
				new_row[1] = tr_start-2 #???
				new_row[3] = ref_seq2
				new_row[4] = alt_seq2

				# vcf_df = vcf_df.append(new_row, ignore_index = True)
				vcf_df.loc[len(vcf_df)] = new_row

			else :
				ref_seq = fa_seq2.seq[tr_start-2]
				ref = str(ref_seq)

				alt = str(ref_seq) + str(alt_seq)
			
			set_ref_alt(ref, alt, i, vcf_df)

	# remove unnecessary columns for VCF
	vcf_df["ID"] = "."

	# adjust variant start position to include reference
	vcf_df["POS"] = vcf_df["POS"] - 1 ### Attention aux translocation avec le nouveau start ?

	# output VCF
	vcf_df.to_csv(output_file, sep="\t", header=False, index=False)



# vcf_header="##source=VISOR BED to VCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"