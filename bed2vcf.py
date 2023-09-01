import pandas as pd
import re, random
from Bio import SeqIO

# read FASTA as a BioSeq dictionary
# (not directly loaded in memory)
def read_fa(input_file):
	record_dict = SeqIO.index(input_file, "fasta")
	return record_dict

def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

def is_reverse(s):
	return(s=="reverse")

def reverse(s):
	return(str(s)[::-1])

def deletion(ref_seq, alt_seq):
	ref = str(ref_seq) + str(alt_seq)
	alt = str(ref_seq)
	return(ref,alt)

def append_li(tup, li1, li2):
	li1.append(tup[0])
	li2.append(tup[1])


# get the sequence of each variants in BED
# output file is a VCF
# (currently, not all variant supported by VISOR are included)
# TODO : missing VISOR variant type : SNP, MNP and 4 types of tandem repeat
def get_seq(vcf_df, bed_df, fa_dict):
	# number of variants
	n = len(vcf_df)
	
	# alternative and reference sequences
	alt = []
	ref = []
	# new variants to add
	rows_to_add = []

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
			ref.append(str(ref_seq) + str(alt_seq))
			alt.append(str(ref_seq))
		# INS
		elif t == "insertion":
			# générer la séquence insérée
			alt_seq = DNA(end-start)
			ref.append(str(ref_seq))
			alt.append(str(ref_seq) + str(alt_seq))
			
		# INV
		elif t == "inversion":
			ref.append(str(ref_seq) + str(alt_seq))
			# reverse alternative sequence
			alt_seq = reverse(alt_seq)
			alt.append(str(ref_seq) + str(alt_seq))

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
			ref.append(str(ref_seq))
			alt.append(str(ref_seq) + str(sv_seq))
			
		# TRA
		# there are 3 translocation types
		elif re.search("translocation", t):
			# get information field
			infos = sv[4]
			# retrieve each info
			info = infos.split(":")
			# translocation start
			tr_start = int(info[2])

			if info[3] == "reverse":
				# reverse alternative sequence
				alt_seq = reverse(alt_seq)
			
			if t == "translocation cut-paste":
				# cut-paste = INS + DEL

				# DEL at the "cut" (reference) position
				ref_seq1 = str(ref_seq) + str(alt_seq)
				alt_seq1 = str(ref_seq)
				ref.append(ref_seq1)
				alt.append(alt_seq1)

				# dupliquer la ligne pour ajouter une insertion
				# new row
				new_sv = sv.tolist()
				# translocation position
				new_sv[1] = tr_start
				# INS at the "paste" position
				ref_seq2 = fa_seq.seq[tr_start-2]
				alt_seq2 = str(ref_seq2) + str(alt_seq)

				new_sv.extend([ref_seq2, alt_seq2])
				rows_to_add.append(new_sv)
				
			elif t == "reciprocal translocation":
				# TODO : gérer orientation --> 2 orientations !!
				# à partir de la position de réf
				# translocation 1
				ref_seq1 = str(ref_seq) + str(alt_seq)
				tr_seq = fa_seq.seq[tr_start-1:tr_start+len(alt_seq)]

				alt_seq1 = str(ref_seq) + str(tr_seq)

				ref.append(ref_seq1)
				alt.append(alt_seq1)
				
				# translocation 2
				# à partir de la position de alt
				ref_tr_seq = fa_seq.seq[tr_start-2]
				ref_seq2 = str(ref_tr_seq) + str(tr_seq)


				if info[4] == "reverse":
					# reverse alternative sequence
					alt_seq = reverse(alt_seq)
	
				alt_seq2 = str(ref_tr_seq) + str(alt_seq)
				# new row
				new_sv = sv.tolist()
				new_sv[1] = tr_start
				new_sv.extend([ref_seq2, alt_seq2])
				rows_to_add.append(new_sv)

			# TODO :  gérer si c'est sur le même chr ou pas = récupérer une autre séquence FASTA dans le dictionnaire
			else :
				ref_seq = fa_seq.seq[tr_start-2]
				ref.append(str(ref_seq))

				alt.append(str(ref_seq) + str(alt_seq))
			
	# create pandas df with reference and alternative sequences
	df = bed_df.assign(ref=ref, alt=alt)
	nrows = pd.DataFrame(rows_to_add, columns=df.columns)
	df = pd.concat([df, nrows], ignore_index=True)
	# remove unnecessary columns for VCF
	df = df.drop(columns=["end","type","info",	"breakpoint"])
	# insert ID column
	n2 = len(df)
	df.insert(2, "id", ["."]*n2)
	# insert other VCF columns with default
	df["qual"] = 99
	df["filter"] = "."
	df["info"] = "."
	df["format"] = "GT"
	df["sample"] = "./."
	# adjust variant start position to include reference
	df["start"] = df["start"] - 1

	# print(df)

	# output VCF
	df.to_csv("out.vcf", sep="\t", header=False, index=False)



vcf_header="##source=VISOR BED to VCF\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"