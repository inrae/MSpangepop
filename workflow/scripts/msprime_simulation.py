# -*- coding: utf-8 -*-
"""
Author: Sukanya Denni, Lucien Piat
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak

This script generates a JSON file for a specified chromosome using msprime simulations.
It utilizes recombination and mutation rates provided by the user and outputs the result to a specified directory.

Usage:
    python script_name.py -fai reference.fai -p 1000 -m 1e-8 -r 1e-8 -n 10 -o output_directory -c chromosome_name
"""

import pandas as pd
import msprime
import argparse
import os
import json
import time
import sys

def get_chromosome_bounds(chrom_length):
    """Define the chromosome boundaries based on length."""
    return [0, chrom_length]

def create_recombination_map(chrom_length, recombination_rate):
    """Create a recombination rate map for the chromosome."""
    try:
        chrom_positions = get_chromosome_bounds(chrom_length)
        return msprime.RateMap(position=chrom_positions, rate=[recombination_rate])
    except Exception as e:
        print(f"MSpangepop -> Error creating recombination map: {e}", file=sys.stderr)
        sys.exit(1)

def save_output(ts_chrom, chromosome_name, output_dir="results"):
    """
    Save simulation results into a JSON file.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)

        # Initialize mutation linkage dictionary
        tree_mutations = {tree_index: [] for tree_index, _ in enumerate(ts_chrom.trees())}

        for mutation in ts_chrom.mutations():
            site = ts_chrom.site(mutation.site)
            mutation_data = {
                "site": mutation.site,
                "site_position": site.position,
                "node": mutation.node,
                "time": mutation.time,
            }
            for tree_index, tree in enumerate(ts_chrom.trees()):
                if tree.interval[0] <= site.position < tree.interval[1]:
                    tree_mutations[tree_index].append(mutation_data)
                    break

        trees_data = []
        for tree_index, tree in enumerate(ts_chrom.trees()):
            edges = [{"parent": parent, "child": child} for parent in tree.nodes() for child in tree.children(parent)]
            tree_data = {
                "tree_index": tree_index,
                "interval": tree.interval,
                "nodes": [{"id": node, "time": ts_chrom.node(node).time} for node in tree.nodes()],
                "edges": edges,
                "mutations": tree_mutations[tree_index],
            }
            trees_data.append(tree_data)

        json_filename = os.path.join(output_dir, f"chr_{chromosome_name}_msprime_simulation.json")
        with open(json_filename, 'w') as f:
            json.dump(trees_data, f, indent=4)

    except Exception as e:
        print(f"MSpangepop -> Error saving output: {e}", file=sys.stderr)
        sys.exit(1)

def simulate_chromosome_evolution(fai_file, population_size, mutation_rate, recombination_rate, sample_size, output_dir, chromosome_name):
    """
    Main function to generate a .json for a specified chromosome using msprime simulations.
    """
    try:
        start_time = time.time()
        
        if not os.path.exists(fai_file):
            raise FileNotFoundError(f"FAI file not found: {fai_file}")

        print("MSpangepop -> Gathering length for chromosome:", chromosome_name)
        chrom_lengths = pd.read_table(fai_file, header=None, usecols=[1], names=["length"])["length"].values
        
        if not chromosome_name.isdigit() or int(chromosome_name) - 1 >= len(chrom_lengths):
            raise ValueError(f"MSpangepop -> Invalid chromosome number: {chromosome_name}")
        
        chrom_length = chrom_lengths[int(chromosome_name) - 1]
        print(f"MSpangepop -> Found chr {chromosome_name} of length: {chrom_length}, starting MSprime simulation...")

        recombination_map = create_recombination_map(chrom_length, recombination_rate)

        ancestry_ts = msprime.sim_ancestry(
            samples=sample_size,
            recombination_rate=recombination_map,
            population_size=population_size
        ).simplify()
        
        mutated_ts = msprime.sim_mutations(ancestry_ts, rate=mutation_rate, discrete_genome=True)

        ts_chrom = mutated_ts.keep_intervals([[0, chrom_length]], simplify=False).trim()

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"MSpangepop -> Simulation for chromosome {chromosome_name} completed in {elapsed_time/60:.2f} min.")

        print(ts_chrom.draw_text())

        os.makedirs(output_dir, exist_ok=True)
        print(f"MSpangepop -> Saving output for chromosome {chromosome_name} to {output_dir}/")
        plt = ts_chrom.draw_svg()
        with open(os.path.join(output_dir, f"chr_{chromosome_name}_tree.svg"), "w") as f:
            f.write(plt)

        save_output(ts_chrom, chromosome_name, output_dir)

    except FileNotFoundError as e:
        print(f"MSpangepop -> File error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"MSpangepop -> Value error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"MSpangepop -> Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a .json for a specific chromosome using msprime simulations.")
    parser.add_argument('-fai', '--fai', type=str, required=True, help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-p', '--population_size', type=int, required=True, help='Effective population size (Ne).')
    parser.add_argument('-m', '--mutation_rate', type=float, required=True, help='Mutation rate per base pair (µ).')
    parser.add_argument('-r', '--recombination_rate', type=float, required=True, help='Recombination rate per base pair.')
    parser.add_argument('-n', '--sample_size', type=int, required=True, help='Sample size (number of individuals to simulate).')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output file.')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='Chromosome to simulate.')

    args = parser.parse_args()

    simulate_chromosome_evolution(
        fai_file=args.fai,
        population_size=args.population_size,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        sample_size=args.sample_size,
        output_dir=args.output_dir,
        chromosome_name=args.chromosome
    )
