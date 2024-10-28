# -*- coding: utf-8 -*-
"""
Author: Sukanya Denni, Lucien Piat
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak

This script generates a VCF file for a specified chromosome using msprime simulations.
It utilizes recombination and mutation rates provided by the user and outputs the result to a specified directory.

Usage:
    python script_name.py -fai reference.fai -p 1000 -m 1e-8 -r 1e-8 -n 10 -o output_directory -c chromosome_name
"""

import pandas as pd
import msprime
import argparse
import os

def get_chromosome_bounds(chrom_length):
    """
    Define the chromosome boundaries (start and end positions) based on length.

    :param chrom_length: int, the length of the chromosome.
    :return: list, [start, end] positions for a single chromosome.
    """
    return [0, chrom_length]

def create_recombination_map(chrom_length, recombination_rate):
    """
    Create a recombination rate map for the chromosome.

    :param chrom_length: int, length of the chromosome.
    :param recombination_rate: float, recombination rate per base pair.
    :return: msprime.RateMap, a map object defining recombination rates across the chromosome.
    """
    chrom_positions = get_chromosome_bounds(chrom_length)
    return msprime.RateMap(position=chrom_positions, rate=[recombination_rate])

def save_vcf_output(ts_chrom, chromosome_name, output_dir="results"):
    """
    Save the VCF file generated from the simulation to the specified directory.

    :param ts_chrom: msprime.TreeSequence, simulated tree sequence with mutations.
    :param chromosome_name: str, the name of the chromosome.
    :param output_dir: str, directory where the VCF file will be saved.
    """
    os.makedirs(output_dir, exist_ok=True)

    vcf_filename = os.path.join(output_dir, f"{chromosome_name}_msprime_simulation.vcf")
    with open(vcf_filename, "w") as vcf_file:
        vcf_file.write(ts_chrom.as_vcf(contig_id=chromosome_name))

def simulate_chromosome_vcf(fai_file, population_size, mutation_rate, recombination_rate, sample_size, output_dir, chromosome_name):
    """
    Main function to generate a VCF for a specified chromosome using msprime simulations.

    :param fai_file: str, path to the FAI index file of the reference FASTA.
    :param population_size: int, effective population size for the simulation.
    :param mutation_rate: float, mutation rate per base pair.
    :param recombination_rate: float, recombination rate per base pair.
    :param sample_size: int, number of samples to simulate.
    :param output_dir: str, directory to save the output VCF file.
    :param chromosome_name: str, specific chromosome to simulate.
    """

    fai_df = pd.read_table(fai_file, header=None, usecols=[0, 1], names=["chromosome", "length"])

    chrom_length = fai_df['length'].values[0]

    recombination_map = create_recombination_map(chrom_length, recombination_rate)

    ancestry_ts = msprime.sim_ancestry(samples=sample_size, recombination_rate=recombination_map, population_size=population_size)
    mutated_ts = msprime.sim_mutations(ancestry_ts, rate=mutation_rate, discrete_genome=True)

    ts_chrom = mutated_ts.keep_intervals([[0, chrom_length]], simplify=False).trim()
    save_vcf_output(ts_chrom, chromosome_name, output_dir)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Generate a VCF for a specific chromosome using msprime simulations.")
    parser.add_argument('-fai', '--fai', type=str, required=True, help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-p', '--population_size', type=int, required=True, help='Effective population size (Ne).')
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help='Mutation rate per base pair (µ).')
    parser.add_argument('-r', '--recombination_rate', type=float, required=True, help='Recombination rate per base pair.')
    parser.add_argument('-n', '--sample_size', type=int, required=True, help='Sample size (number of individuals to simulate).')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output VCF file.')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='Chromosome to simulate.')

    args = parser.parse_args()

    simulate_chromosome_vcf(
        fai_file=args.fai,
        population_size=args.population_size,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        sample_size=args.sample_size,
        output_dir=args.output_dir,
        chromosome_name=args.chromosome
    )
