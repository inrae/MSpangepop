# fusionner le VCF de msprime (position, génotype) et le BED de VISOR
# pour obtenir un BED avec tous les variants
# récupérer les variants et les génotypes msprime
import io, os, shutil
import pandas as pd

# def read_vcf(input_file, header_size=6):
# 	with open(input_file, "r") as txt:
# 		t = txt.read()
# 	txt.close()

# 	# write header for later
# 	if not os.path.exists("results"):
# 		os.mkdir("results")
# 	f = open("results/vcf_header.txt", "w")
# 	f.write(''.join(t.splitlines(keepends=True)[:header_size]))
# 	f.close()

# 	# remove VCF header for dataframe
# 	t = ''.join(t.splitlines(keepends=True)[header_size-1:])
# 	df = pd.read_table(io.StringIO(t))
# 	df = df.rename(columns={"#CHROM" : "CHROM"})
# 	return(df)

def get_vcf(input_file):
	df = pd.read_table(input_file, sep="\t")
	return(df)

# def read_BED(input_file):
# 	colnames=["chr", "start", "end", "type","info","breakpoint"]
# 	df = pd.read_table(input_file, header=None, names=colnames)
# 	return(df)

# def replace_bed_col(input_BED, input_VCF):
# 	vcf = read_vcf(input_VCF)
# 	bed = read_BED(input_BED)
# 	if len(bed) == len(vcf):
# 		bed["chr"] = vcf["CHROM"]
# 		bed["start"] = vcf["POS"]
# 	return(bed)

# def merge_final_vcf(header_file, content_file, merged_file):
#     with open(merged_file, 'wb') as outfile:
#         for filename in [header_file, content_file]:
#             with open(filename, 'rb') as infile:
#                 shutil.copyfileobj(infile, outfile)