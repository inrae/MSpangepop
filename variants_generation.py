from bed2vcf import read_fa, get_seq
from fusion import *
import argparse
from multiprocessing import Process

## n : length of variant
def sv_vcf(bed, vcf, fasta, outName, header):
	bed_df = replace_bed_col(bed, vcf)
	if header != None:
		vcf_df=read_vcf(vcf, header)
	else:
		vcf_df=read_vcf(vcf)
	fa = read_fa(fasta)

	get_seq(vcf_df, bed_df, fa, outName)

	header_file="results/vcf_header.txt"
	content_file = "results/" + outName + ".vcf"
	merged_file="results/" + outName + "_final.vcf"
	merge_final_vcf(header_file, content_file, merged_file)

# parse arguments
parser = argparse.ArgumentParser(description='create final VCF')
parser.add_argument('-b','--bed', type=str, required = True, help='BED file of structural variants to include')
parser.add_argument('-v', '--vcf', type=str, required = True, help='VCF generated with msprime --> tree generation script')
parser.add_argument('-fa', '--fasta', type=str, required = True, help='Reference FASTA for the VCF and the variants')
parser.add_argument('-o', '--outName', type=str, required = True, help='Output name')
parser.add_argument('--header', type=int, help='vcf header size')

if __name__ == '__main__':
	args = parser.parse_args()
	d = vars(args)
	li = list(d.values())
	p = Process(target=sv_vcf, args=(li))
	p.start()
	p.join()