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
        print(f"❌ MSpangepop -> Error creating recombination map: {e}", file=sys.stderr)
        sys.exit(1)

def save_output(mutated_ts, chromosome_name, output_dir="results"):
    """
    Save simulation results into a JSON file efficiently for large files.
    """
    try:
        os.makedirs(output_dir, exist_ok=True)

        json_filename = os.path.join(output_dir, f"chr_{chromosome_name}_msprime_simulation.json")

        with open(json_filename, 'w') as f:
            f.write('[\n')  # Start JSON array

            first_tree = True
            for tree_index, tree in enumerate(mutated_ts.trees()):
                # Convert interval to integers
                interval = [int(tree.interval[0]), int(tree.interval[1])]

                edges = [{"parent": parent, "child": child} for parent in tree.nodes() for child in tree.children(parent)]
                nodes = [{"id": node, "time": mutated_ts.node(node).time} for node in tree.nodes()]

                mutations = []
                for mutation in mutated_ts.mutations():
                    site = mutated_ts.site(mutation.site)
                    if tree.interval[0] < site.position < tree.interval[1]:
                        mutations.append({
                            "site_position": int(site.position),  # Convert to int
                            "node": mutation.node,
                            "time": mutation.time
                        })
                tree_data = {
                    "tree_index": tree_index,
                    "interval": interval,
                    "nodes": nodes,
                    "edges": edges,
                    "mutations": mutations,
                }

                if not first_tree:
                    f.write(",\n")  # Separate JSON objects
                first_tree = False

                json.dump(tree_data, f, indent=4)

            f.write('\n]')  # End JSON array

    except Exception as e:
        print(f"❌ MSpangepop -> Error saving output: {e}", file=sys.stderr)
        sys.exit(1)


def simulate_chromosome_evolution(fai_file, population_size, mutation_rate, recombination_rate, sample_size, output_dir, chromosome_name, model):
    """
    Main function to generate a .json for a specified chromosome using msprime simulations.
    """
    try:
        start_time = time.time()
        
        if not os.path.exists(fai_file):
            raise FileNotFoundError(f"❌ MSpangepop -> FAI file not found: {fai_file}")

        print("🔹 MSpangepop -> Gathering length for chromosome:", chromosome_name)
        chrom_lengths = pd.read_table(fai_file, header=None, usecols=[1], names=["length"])["length"].values
        
        if not chromosome_name.isdigit() or int(chromosome_name) - 1 >= len(chrom_lengths):
            raise ValueError(f"❌ MSpangepop ->  Invalid chromosome number: {chromosome_name}")
        
        chrom_length = chrom_lengths[int(chromosome_name) - 1]
        print(f"✅ MSpangepop -> Found chr {chromosome_name} of length: {chrom_length},\n\t🔹 Starting MSprime simulation...")

        recombination_map = create_recombination_map(chrom_length, recombination_rate)

        ancestry_ts = msprime.sim_ancestry(
            samples=sample_size,
            recombination_rate=recombination_map,
            population_size=population_size
        ).simplify()

        plt = ancestry_ts.draw_svg()
        with open(os.path.join(output_dir, f"chr_{chromosome_name}_ancestery.svg"), "w") as f:
            f.write(plt)

        mutated_ts = msprime.sim_mutations(ancestry_ts, rate=mutation_rate, discrete_genome=True, model=model)
        mutated_ts = mutated_ts.keep_intervals([[0, chrom_length]], simplify=True).trim()


        print(f"✅ MSpangepop -> Simulation for chromosome {chromosome_name}")

        print(mutated_ts.draw_text())

        os.makedirs(output_dir, exist_ok=True)
        print(f"🔹 MSpangepop -> Saving output for chromosome {chromosome_name} to {output_dir}/")
        
        plt = mutated_ts.draw_svg()
        with open(os.path.join(output_dir, f"chr_{chromosome_name}_mutations.svg"), "w") as f:
            f.write(plt)

        save_output(mutated_ts, chromosome_name, output_dir)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"✅ MSpangepop -> Simulation saved for chromosome {chromosome_name} completed in {elapsed_time/60:.2f} min.")
        

    except FileNotFoundError as e:
        print(f"❌ MSpangepop -> File error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"❌ MSpangepop -> Value error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ MSpangepop -> Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a .json for a specific chromosome using msprime simulations.")
    parser.add_argument('-fai', '--fai', type=str, required=True, help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-p', '--population_size', type=int, required=True, help='Effective population size (Ne).')
    parser.add_argument('-mu', '--mutation_rate', type=float, required=True, help='Mutation rate per base pair (µ).')
    parser.add_argument('-r', '--recombination_rate', type=float, required=True, help='Recombination rate per base pair.')
    parser.add_argument('-n', '--sample_size', type=int, required=True, help='Sample size (number of individuals to simulate).')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output file.')
    parser.add_argument('-c', '--chromosome', type=str, required=True, help='Chromosome number (1-based index).')
    parser.add_argument('-mo', '--model', type=str, required=True, help='msprime mutation model')
    args = parser.parse_args()

    simulate_chromosome_evolution(
        fai_file=args.fai,
        population_size=args.population_size,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        sample_size=args.sample_size,
        output_dir=args.output_dir,
        chromosome_name=args.chromosome,
        model=args.model
    )
