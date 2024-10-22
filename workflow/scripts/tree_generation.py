# -*- coding: utf-8 -*-
import pandas as pd
import math
import msprime
import argparse
import os

def chrom_pos(length):
    """Calculate chromosome positions based on length."""
    return [0, length]  # For a single chromosome, just return 0 and its length

def chrom_rate():
    """Return a fixed rate map for the single chromosome."""
    r_chrom = 1e-8  # Recombination rate for the entire chromosome
    return [r_chrom]  # Only one rate for the single interval

def msprime_map(length):
    """Create a rate map for msprime based on chromosome length."""
    map_positions = chrom_pos(length)
    rates = chrom_rate()
    return msprime.RateMap(position=map_positions, rate=rates)

def output_files(ts_chrom, name, outdir="results"):
    """Output VCF files based on the simulation results."""
    os.makedirs(outdir, exist_ok=True)  # Create output directory if it doesn't exist

    # Write VCF for the single chromosome
    filename = os.path.join(outdir, f"{name}_msprime_simulation.vcf")
    with open(filename, "w") as vcf_file:
        vcf_file.write(ts_chrom.as_vcf(contig_id=name))

def msprime_vcf(fai, pop_size, mut_rate, n, outdir, chromosome):
    """Main function to generate VCF for a specific chromosome using msprime."""
    df = pd.read_table(fai, header=None, usecols=[0, 1], names=["name", "length"])
    
    # Filter for the specific chromosome
    chrom_data = df[df['name'] == chromosome]
    if chrom_data.empty:
        raise ValueError(f"Chromosome {chromosome} not found in FAI file.")

    length = chrom_data['length'].values[0]  # Get the length of the specific chromosome
    rate_map = msprime_map(length)

    # Simulate ancestry and mutations
    ts = msprime.sim_ancestry(samples=n, recombination_rate=rate_map, population_size=pop_size)
    mutated_ts = msprime.sim_mutations(ts, rate=mut_rate, discrete_genome=True)

    # Keep the interval for the specific chromosome
    ts_chrom = mutated_ts.keep_intervals([[0, length]], simplify=False).trim()
    output_files(ts_chrom, chromosome, outdir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate VCF for a specific chromosome in reference FASTA.')
    parser.add_argument('-fai', '--fai', type=str, required=True, help='FAI samtools index of reference FASTA')
    parser.add_argument('-p', '--popSize', type=int, required=True, help='population size (Ne)')
    parser.add_argument('-r', '--rate', type=float, required=True, help='mutation rate (µ)')
    parser.add_argument('-n', '--sampleSize', type=int, required=True, help='sample size')
    parser.add_argument('-o', '--outDir', type=str, required=True, help='output directory')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='specific chromosome to simulate')

    args = parser.parse_args()

    msprime_vcf(args.fai, args.popSize, args.rate, args.sampleSize, args.outDir, args.chromosome)
