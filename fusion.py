# fusionner le VCF de msprime (position, génotype) et le BED de VISOR
# pour obtenir un BED avec tous les variants
# récupérer les variants et les génotypes msprime
import io
import pandas as pd

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
		bed["end"] = vcf["POS"] + len_SV
	return(bed)

# étape 1 : vérifier que les 2 tableaux font la même taille
# 2 : prendre les positions du VCF + taille = end
# 3 : remplacer les colonnes start end par les positions VCF
# 4 : enlever le breakpoint ?
# 5 avec le fasta, récupérer les séquences de variants
# 6 : créer le nouveau VCF -> reprendre le VCF msprime et coller les bon variants à la place