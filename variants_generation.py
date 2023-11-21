from bed2vcf import read_fa, get_seq
from randombed import generate_type
import argparse
import pandas as pd
from multiprocessing import Process

## n : length of variant
def sv_vcf(vcf, fasta, fai, yml, outName):
	vcf_df = pd.read_table(vcf, sep="\t")
	nvar = len(vcf_df)
	bed_df = generate_type(nvar, yml, fai)
	fa = read_fa(fasta)

	get_seq(vcf_df, bed_df, fa, outName)

# parse arguments
parser = argparse.ArgumentParser(description='create final VCF')
# parser.add_argument('-b','--bed', type=str, required = True, help='BED file of structural variants to include')
parser.add_argument('-v', '--vcf', type=str, required = True, help='VCF generated with msprime --> tree generation script')
parser.add_argument('-fa', '--fasta', type=str, required = True, help='Reference FASTA for the VCF and the variants')
parser.add_argument('-fai', '--fai', type=str, required = True, help='Samtools index of reference FASTA')
parser.add_argument('-y', '--yaml', type=str, required = True, help='Reference FASTA for the VCF and the variants')
parser.add_argument('-o', '--outName', type=str, required = True, help='Output name')
# parser.add_argument('--header', type=int, help='vcf header size')

if __name__ == '__main__':
	args = parser.parse_args()
	d = vars(args)
	li = list(d.values())
	p = Process(target=sv_vcf, args=(li))
	p.start()
	p.join()