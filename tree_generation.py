import pandas as pd
import math, msprime, argparse, os
import numpy as np
from multiprocessing import Process

def chrom_pos(df):
    l=df["length"].to_list()
    map_positions=[0]
    for i in l:
        map_positions.append(i+map_positions[-1])
    n=len(map_positions)
    
    for i in range(1,n*2-3, 2):
        map_positions.insert(i+1, map_positions[i]+1)
    return(map_positions)


### simulation de plusieurs chromosomes
def chrom_rate(df):
    n=len(df)
    r_chrom = 1e-8
    r_break = math.log(2)
    li = np.resize([r_chrom,r_break], (n*2-1))
    return(li)

def msprime_map(df):
    map_positions=chrom_pos(df)
    rates=chrom_rate(df)
    rate_map=msprime.RateMap(position=map_positions, rate=rates)
    return(rate_map)

### ts_chroms = list of trees for each chromosome
### names = list of chromosome name
def output_files(ts_chroms, names, outdir="results"):
    if not os.path.exists(outdir):
            os.mkdir(outdir)

        ## write vcf
    filename = outdir + "/msprime_simulation.vcf"
    for i in range(len(ts_chroms)):
        txt = ts_chroms[i].as_vcf(contig_id=names[i])
        t = ''.join(txt.splitlines(keepends=True)[6:])
        # filename = 'tmp_vcf/output' + str(i) + ".vcf"
        if i==0:
            header_top = ''.join(txt.splitlines(keepends=True)[:4])
            header_bot = ''.join(txt.splitlines(keepends=True)[4:6])
            with open(filename, "w") as vcf_file:
                vcf_file.write(t)
                # ts_chroms[i].write_vcf(vcf_file, contig_id=names[i])
        else:
            header_top += txt.splitlines(keepends=True)[3]
            
            with open(filename, "a") as vcf_file:
                vcf_file.write(t)
            vcf_file.close()
    ## write header to file
    header = header_top + header_bot
    with open(outdir + "/header.vcf", "w") as head:
        head.write(header)
    head.close()

# input_fai = samtools FAI of FASTA where variants will be generated
# pop_size = population size
# mut_rate = mutation rate
# n = sample size / number of indiv
def msprime_vcf(fai, pop_size, mut_rate, n, outdir):

    df = pd.read_table(fai, header=None, usecols=[0,1], names =["name", "length"])
    rate_map=msprime_map(df)

    ts = msprime.sim_ancestry(
		samples=n,
		recombination_rate=rate_map,
		population_size=pop_size,
		# random_seed=123456
        )

	# avec des mutations
    mutated_ts = msprime.sim_mutations(
        ts, rate=mut_rate, 
        # random_seed=5678, 
        discrete_genome=True)

    ts_chroms = []

    l=df["length"].to_list()

    chrom_positions=[0]
    for i in l:
        chrom_positions.append(i+chrom_positions[-1])

    for j in range(len(chrom_positions) - 1):
        start, end = chrom_positions[j: j + 2]
        chrom_ts = mutated_ts.keep_intervals([[start, end]], simplify=False).trim()
        ts_chroms.append(chrom_ts)

    names=df["name"].to_list()

    ## create files    
    if outdir == None:
        output_files(ts_chroms, names)
    else:
        output_files(ts_chroms, names, outdir)

parser = argparse.ArgumentParser(description='Generate VCF for each chromosome in reference FASTA.')
parser.add_argument('-fai', '--fai', type=str, required = True, help='FAI samtools index of reference FASTA')
parser.add_argument('-p', '--popSize', type=int, required = True, help='population size (Ne)')
parser.add_argument('-r', '--rate', type=float, required = True, help='mutation rate (µ)')
parser.add_argument('-n', '--sampleSize', type=int, required = True, help='sample size')
parser.add_argument('-o', '--outDir', type=str, help='output directory')

if __name__ == '__main__':
	args = parser.parse_args()
	d = vars(args)
	li = list(d.values())
	p = Process(target=msprime_vcf, args=(li))
	p.start()
	p.join()