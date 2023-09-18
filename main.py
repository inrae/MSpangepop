from bed2vcf import read_fa, get_seq
from fusion import *
from tree_generation import msprime_vcf

#n : length of variant
def main(bed, vcf, fasta, n, output_file):
	
	bed_df = replace_bed_col(bed, VCF, n)
	vcf_df=read_vcf(vcf)
	fa = read_fa(fasta)


	get_seq(vcf_df, bed_df, fa, output_file)

bed = "mini_random.bed"
vcf = "output.vcf"
fasta = "/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/g1.chr8.fasta"
n = 200

output_file="60_tritici.vcf"


# generate msprime population VCF
fai = "/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/g1.chr8.fasta.fai"

# samtools FAI
# population size
# mutation rate
msprime_vcf(fai, 1000, 1e-8)