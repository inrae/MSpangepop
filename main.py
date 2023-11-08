from bed2vcf import read_fa, get_seq
from fusion import *
from tree_generation import msprime_vcf
import defopt

## n : length of variant
def sv_vcf(bed: str, vcf: str, fasta: str, outName: str):
	"""
    Create a VCF with structural variants.
    
    :parameter bed: BED file of structural variants to include
    :parameter vcf: VCF generated with msprime --> tree generation script
	:parameter fasta: Reference FASTA for the VCF and the variants
	:parameter outName: Output name
    
    """
	bed_df = replace_bed_col(bed, vcf)
	vcf_df=read_vcf(vcf)
	fa = read_fa(fasta)

	get_seq(vcf_df, bed_df, fa, outName)

	header_file="results/vcf_header.txt"
	content_file = "results/" + outName + ".vcf"
	merged_file="results/" + outName + "_final.vcf"
	merge_final_vcf(header_file, content_file, merged_file)

if __name__ == '__main__':
    defopt.run([sv_vcf, msprime_vcf])

bed = "test/mini_random.bed"
vcf = "test/output.vcf"
fasta = "/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/ztIPO323.chr8.fasta"

output_file="tritici_test"
# main(bed, vcf, fasta, output_file)

### generate msprime population VCF
# fai = "ztIPO323.chr8.fasta.fai"

# ## samtools FAI
# ## population size
# ## mutation rate
# ## sample size
# msprime_vcf(fai, 5000, 1e-8, 60)