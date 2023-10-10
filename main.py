from bed2vcf import read_fa, get_seq
from fusion import *
from tree_generation import msprime_vcf

## n : length of variant
def main(bed, vcf, fasta, output_file):
	
	bed_df = replace_bed_col(bed, vcf)
	vcf_df=read_vcf(vcf)
	fa = read_fa(fasta)


	get_seq(vcf_df, bed_df, fa, output_file)

bed = "test/mini_random.bed"
vcf = "test/output.vcf"
fasta = "/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/ztIPO323.chr8.fasta"

output_file="tritici_test"
main(bed, vcf, fasta, output_file)

header_file="results/vcf_header.txt"
content_file = "results/" + output_file + ".vcf"
merged_file="results/" + output_file + "_final.vcf"
merge_final_vcf(header_file, content_file, merged_file)

### generate msprime population VCF
# fai = "ztIPO323.chr8.fasta.fai"

# ## samtools FAI
# ## population size
# ## mutation rate
# ## sample size
# msprime_vcf(fai, 5000, 1e-8, 60)