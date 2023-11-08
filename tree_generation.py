import pandas as pd
import math, msprime
import numpy as np

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

# input_fai = samtools FAI of FASTA where variants will be generated
# pop_size = population size
# mut_rate = mutation rate
# n = sample size / number of indiv
def msprime_vcf(fai: str, pop_size: int, mut_rate: float, n: int):
    """
    Generate VCF for each chromosome in reference FASTA.

    :parameter fai: FAI samtools index of reference FASTA
    :parameter pop_size: population size (Ne)
    :parameter mut_rate: mutation rate (µ)
    :parameter n: sample size
    """
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
        discrete_genome=False)

    ts_chroms = []

    l=df["length"].to_list()

    chrom_positions=[0]
    for i in l:
        chrom_positions.append(i+chrom_positions[-1])

    for j in range(len(chrom_positions) - 1):
        start, end = chrom_positions[j: j + 2]
        chrom_ts = mutated_ts.keep_intervals([[start, end]], simplify=False).trim()
        ts_chroms.append(chrom_ts)
        print(chrom_ts.sequence_length)

    names=df["name"].to_list()

    for i in range(len(ts_chroms)):
        filename="msprime_output_chr.vcf"
        if i==0:
            with open(filename, "w") as vcf_file:
                ts_chroms[i].write_vcf(vcf_file, contig_id=names[i])
        else:
            txt = ts_chroms[i].as_vcf(contig_id=names[i])
            t = ''.join(txt.splitlines(keepends=True)[6:])
            with open(filename, "a") as vcf_file:
                vcf_file.write(t)
            vcf_file.close()

    print(len(ts_chroms))