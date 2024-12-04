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
import json

def get_chromosome_bounds(chrom_length):
    """
    Define the chromosome boundaries based on length.
    """
    return [0, chrom_length]

def create_recombination_map(chrom_length, recombination_rate):
    """
    Create a recombination rate map for the chromosome.
    """
    chrom_positions = get_chromosome_bounds(chrom_length)
    return msprime.RateMap(position=chrom_positions, rate=[recombination_rate])
import json

def save_output(ts_chrom, chromosome_name, output_dir="results", batch_size=1000):
    """
    Save all simulation trees and mutation information into a single JSON file, 
    linking mutations to the trees they occur in.

    :param ts_chrom: Tree sequence object
    :param chromosome_name: Name of the chromosome simulated
    :param output_dir: Directory to save the output files
    :param batch_size: Number of records per batch (if needed)
    """
    
    # Initialize mutation linkage dictionary
    tree_mutations = {tree_index: [] for tree_index, _ in enumerate(ts_chrom.trees())}
    
    # Collect mutations and link them to their corresponding trees
    for mutation in ts_chrom.mutations():
        site = ts_chrom.site(mutation.site)
        mutation_data = {
            "site": mutation.site,
            "site_position": site.position,
            "node": mutation.node,
            "time": mutation.time,
        }
        
        # Find the tree containing this mutation
        for tree_index, tree in enumerate(ts_chrom.trees()):
            if tree.interval[0] <= site.position < tree.interval[1]:
                tree_mutations[tree_index].append(mutation_data)
                break
    
    # Collect all tree data
    trees_data = []
    for tree_index, tree in enumerate(ts_chrom.trees()):
        # Retrieve all parent-child relationships in the tree
        edges = []
        for parent in tree.nodes():
            for child in tree.children(parent):
                edges.append({"parent": parent, "child": child})
        
        tree_data = {
            "tree_index": tree_index,
            "interval": tree.interval,
            "nodes": [{"id": node, "time": ts_chrom.node(node).time} for node in tree.nodes()],
            "edges": edges,  # Use computed edges
            "mutations": tree_mutations[tree_index],  # Add mutations specific to this tree
        }
        trees_data.append(tree_data)
    
    # Output as JSON
    json_filename = os.path.join(output_dir, f"{chromosome_name}_msprime_simulation.json")
    with open(json_filename, 'w') as f:
        json.dump(trees_data, f, indent=4)


def simulate_chromosome_vcf(fai_file, population_size, mutation_rate, recombination_rate, sample_size, output_dir, chromosome_name):
    """
    Main function to generate a VCF for a specified chromosome using msprime simulations.
    After simulation, the tree is visualized and saved as a PNG image.

    :param fai_file: str, path to the FAI index file of the reference FASTA.
    :param population_size: int, effective population size for the simulation.
    :param mutation_rate: float, mutation rate per base pair.
    :param recombination_rate: float, recombination rate per base pair.
    :param sample_size: int, number of samples to simulate.
    :param output_dir: str, directory to save the output VCF file.
    :param chromosome_name: str, specific chromosome to simulate.
    """
    # Load the chromosome length from the FAI file (assuming only one chromosome in FAI file)
    chrom_length = pd.read_table(fai_file, header=None, usecols=[1], names=["length"])['length'].values[0]

    recombination_map = create_recombination_map(chrom_length, recombination_rate)

    # Simulate the ancestry and mutations
    ancestry_ts = msprime.sim_ancestry(
        samples=sample_size,
        recombination_rate=recombination_map,
        population_size=population_size
    )
    ancestry_ts = ancestry_ts.simplify()
    mutated_ts = msprime.sim_mutations(ancestry_ts, rate=mutation_rate, discrete_genome=True)
    
    # Trim and keep the relevant intervals
    ts_chrom = mutated_ts.keep_intervals([[0, chrom_length]], simplify=False).trim()
    print(ts_chrom.draw_text())


    plt = ts_chrom.draw_svg()
    with open(f"{output_dir}/{chromosome_name}_tree.svg", "w") as f:
        f.write(plt) # Write the SVG string to a file

    save_output(ts_chrom, chromosome_name, output_dir)


if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate a VCF for a specific chromosome using msprime simulations.")
    parser.add_argument('-fai', '--fai', type=str, required=True, help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-p', '--population_size', type=int, required=True, help='Effective population size (Ne).')
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help='Mutation rate per base pair (µ).')
    parser.add_argument('-r', '--recombination_rate', type=float, required=True, help='Recombination rate per base pair.')
    parser.add_argument('-n', '--sample_size', type=int, required=True, help='Sample size (number of individuals to simulate).')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output VCF file.')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='Chromosome to simulate.')

    args = parser.parse_args()

    # Run the simulation with specified parameters
    simulate_chromosome_vcf(
        fai_file=args.fai,
        population_size=args.population_size,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        sample_size=args.sample_size,
        output_dir=args.output_dir,
        chromosome_name=args.chromosome
    )