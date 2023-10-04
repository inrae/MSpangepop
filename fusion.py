# fusionner le VCF de msprime (position, génotype) et le BED de VISOR
# pour obtenir un BED avec tous les variants
# récupérer les variants et les génotypes msprime
import io, json
import pandas as pd
import numpy as np

def read_vcf(input_file):
	with open(input_file, "r") as txt:
		t = txt.read()
	txt.close()
	t = ''.join(t.splitlines(keepends=True)[5:])
	
	df = pd.read_table(io.StringIO(t))
	df = df.rename(columns={"#CHROM" : "CHROM"})
	return(df)

def read_BED(input_file):
	colnames=["chr", "start", "end", "type","info","breakpoint"]
	df = pd.read_table(input_file, header=None, names=colnames)
	return(df)

def replace_bed_col(input_BED, input_VCF, len_SV):
	vcf = read_vcf(input_VCF)
	bed = read_BED(input_BED)
	if len(bed) == len(vcf):
		bed["chr"] = vcf["CHROM"]
		bed["start"] = vcf["POS"]
		# piocher dans une distribution de taille de SV --> modéliser sur la distrib humain, etc
		bed["end"] = vcf["POS"] + len_SV
	return(bed)

def get_random_len(n):
	df = pd.read_csv("/home/sukanya/tests/03_results/SV_distrib_human/size_distribINS.tsv", sep="\t")
	pb = df["pb"].tolist()
	li = np.random.multinomial(n, pb)
	print(li)
	s = []
	for i in range(len(li)):
	
		if li[i] != 0:
			interval = df["size_interval"].iloc[i]
			interval = json.loads(interval)
			# print(interval)
			s += list(np.random.uniform(interval[0], interval[1], li[i]).round())
	print(s)

get_random_len(10)


# étape 1 : vérifier que les 2 tableaux font la même taille
# 2 : prendre les positions du VCF + taille = end
# 3 : remplacer les colonnes start end par les positions VCF
# 4 : enlever le breakpoint ?
# 5 avec le fasta, récupérer les séquences de variants
# 6 : créer le nouveau VCF -> reprendre le VCF msprime et coller les bon variants à la place