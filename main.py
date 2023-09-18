from bed2vcf import read_fa, get_seq
from fusion import *
from tree_generation import msprime_vcf

## n : length of variant
def main(bed, vcf, fasta, n, output_file):
	
	bed_df = replace_bed_col(bed, vcf, n)
	vcf_df=read_vcf(vcf)
	fa = read_fa(fasta)


	get_seq(vcf_df, bed_df, fa, output_file)

# bed = "mini_random.bed"
# vcf = "output.vcf"
# fasta = "ztIPO323.chr8.fasta"
# n = 200

# output_file="60_tritici.vcf"


### generate msprime population VCF
fai = "ztIPO323.chr8.fasta.fai"

## samtools FAI
## population size
## mutation rate
## sample size
msprime_vcf(fai, 5000, 1e-8, 60)