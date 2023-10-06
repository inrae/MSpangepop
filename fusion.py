# fusionner le VCF de msprime (position, génotype) et le BED de VISOR
# pour obtenir un BED avec tous les variants
# récupérer les variants et les génotypes msprime
import io
import pandas as pd

def read_vcf(input_file):
	with open(input_file, "r") as txt:
		t = txt.read()
	txt.close()
	# remove VCF header for dataframe
	t = ''.join(t.splitlines(keepends=True)[5:])
	df = pd.read_table(io.StringIO(t))
	df = df.rename(columns={"#CHROM" : "CHROM"})
	return(df)

def read_BED(input_file):
	colnames=["chr", "start", "end", "type","info","breakpoint"]
	df = pd.read_table(input_file, header=None, names=colnames)
	return(df)

def replace_bed_col(input_BED, input_VCF):
	vcf = read_vcf(input_VCF)
	bed = read_BED(input_BED)
	if len(bed) == len(vcf):
		bed["chr"] = vcf["CHROM"]
		bed["start"] = vcf["POS"]
		# piocher dans une distribution de taille de SV --> modéliser sur la distrib humain, etc
		# bed["end"] = vcf["POS"] + len_SV
	return(bed)