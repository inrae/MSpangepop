# -*- coding: utf-8 -*-
import pandas as pd
import math
import msprime
import argparse
import os
import numpy as np
from multiprocessing import Process

def chrom_pos(lengths):
    """Calculate chromosome positions based on lengths."""
    map_positions = [0]
    for length in lengths:
        map_positions.append(length + map_positions[-1])
    
    # Insert midpoints for recombination
    for i in range(1, len(map_positions) * 2 - 3, 2):
        map_positions.insert(i + 1, map_positions[i] + 1)
    
    return map_positions

def chrom_rate(num_chromosomes):
    """Return a rate map for chromosomes."""
    r_chrom = 1e-8
    r_break = math.log(2)
    return np.resize([r_chrom, r_break], (num_chromosomes * 2 - 1))

def msprime_map(lengths):
    """Create a rate map for msprime."""
    map_positions = chrom_pos(lengths)
    rates = chrom_rate(len(lengths))
    return msprime.RateMap(position=map_positions, rate=rates)

def output_files(ts_chroms, names, outdir="results"):
    """Output VCF files and headers based on the simulation results."""
    os.makedirs(outdir, exist_ok=True)  # Create output directory if it doesn't exist

    # Write VCF
    filename = os.path.join(outdir, "msprime_simulation.vcf")
    header_top, header_bot = "", ""
    
    for i, ts in enumerate(ts_chroms):
        txt = ts.as_vcf(contig_id=names[i])
        vcf_content = ''.join(txt.splitlines(keepends=True)[6:])  # Extract VCF records

        # Capture header information
        if i == 0:
            header_top = ''.join(txt.splitlines(keepends=True)[:4])
            header_bot = ''.join(txt.splitlines(keepends=True)[4:6])
            with open(filename, "w") as vcf_file:
                vcf_file.write(vcf_content)
        else:
            header_top += txt.splitlines(keepends=True)[3]
            with open(filename, "a") as vcf_file:
                vcf_file.write(vcf_content)

    # Write header to file
    header = header_top + header_bot
    with open(os.path.join(outdir, "header.vcf"), "w") as head_file:
        head_file.write(header)

def msprime_vcf(fai, pop_size, mut_rate, n, outdir):
    """Main function to generate VCFs using msprime."""
    df = pd.read_table(fai, header=None, usecols=[0, 1], names=["name", "length"])
    rate_map = msprime_map(df["length"].to_list())

    ts = msprime.sim_ancestry(samples=n, recombination_rate=rate_map, population_size=pop_size)
    mutated_ts = msprime.sim_mutations(ts, rate=mut_rate, discrete_genome=True)

    ts_chroms = []
    chrom_positions = np.cumsum(df["length"].to_list()).tolist()

    for start, end in zip(chrom_positions[:-1], chrom_positions[1:]):
        chrom_ts = mutated_ts.keep_intervals([[start, end]], simplify=False).trim()
        ts_chroms.append(chrom_ts)

    names = df["name"].to_list()
    output_files(ts_chroms, names, outdir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate VCF for each chromosome in reference FASTA.')
    parser.add_argument('-fai', '--fai', type=str, required=True, help='FAI samtools index of reference FASTA')
    parser.add_argument('-p', '--popSize', type=int, required=True, help='population size (Ne)')
    parser.add_argument('-r', '--rate', type=float, required=True, help='mutation rate (µ)')
    parser.add_argument('-n', '--sampleSize', type=int, required=True, help='sample size')
    parser.add_argument('-o', '--outDir', type=str, help='output directory')

    args = parser.parse_args()
    # Pass arguments as a tuple
    p = Process(target=msprime_vcf, args=(args.fai, args.popSize, args.rate, args.sampleSize, args.outDir))
    p.start()
    p.join()
