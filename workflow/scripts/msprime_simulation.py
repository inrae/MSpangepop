# -*- coding: utf-8 -*-
"""
Author: Sukanya Denni, Lucien Piat
Institution: INRAe
Project: PangenOak

This script generates a JSON file for a specified chromosome using msprime simulations.
It utilizes recombination and mutation rates provided by the user and outputs the result to a specified directory.
"""

import msprime
import argparse
import os
import json
import time

from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, get_indent, process_seed

def get_chromosome_bounds(chrom_length: int) -> tuple:
    """Define the chromosome boundaries based on length."""
    if chrom_length < 0:
        raise MSerror("Chromosome length cannot be negative.")
    return (0, chrom_length)

def create_recombination_map(chrom_length:int, recombination_rate:float)-> msprime.RateMap:
    """Create a recombination rate map for the chromosome."""
    try:
        return msprime.RateMap(position=get_chromosome_bounds(chrom_length), rate=[recombination_rate])
    except Exception as e:
            raise MSerror(f"Error creating recombination map: {e}")

def save_output(mutated_ts, chromosome_name: str, output_dir: str, readable_json: bool):
    """
    Save simulation results into a JSON file efficiently for large files.
    Uses incremental JSON writing to avoid high memory usage.
    """

    json_filename = os.path.join(output_dir, f"chr_{chromosome_name}_msprime_simulation.json")

    try:
        with open(json_filename, 'w') as f:
            f.write('[\n')  # Start JSON array
            first_tree = True

            for tree_index, tree in enumerate(mutated_ts.trees()):
                interval = [int(tree.interval[0]), int(tree.interval[1])]  # Tree interval

                # Cache nodes and children to avoid redundant function calls
                nodes = [{"id": node, "time": mutated_ts.node(node).time} for node in tree.nodes()]
                edges = []
                for parent in tree.nodes():
                    children = list(tree.children(parent))
                    edges.extend({"parent": parent, "child": child} for child in children)

                # Collect mutations within tree interval
                mutations = [
                    {"node": mutation.node}
                    for mutation in mutated_ts.mutations()
                    if tree.interval[0] < mutated_ts.site(mutation.site).position < tree.interval[1]
                ]

                # Construct tree data
                tree_data = {
                    "tree_index": tree_index,
                    "interval": interval,
                    "nodes": nodes,
                    "edges": edges,
                    "mutations": mutations,
                }

                # Write JSON object, separating with a comma if not the first entry
                if not first_tree:
                    f.write(",\n")
                first_tree = False

                json.dump(tree_data, f, indent=get_indent(readable_json))

            f.write('\n]')  # End JSON array

        MSsuccess(f"Successfully saved output to {json_filename}")

    except Exception as e:
        raise MSerror(f"Error saving JSON output: {e}")

def get_chromosome_length(fai_file, chromosome_name):
    """
    Retrieves the length of a specified chromosome from a FASTA index (.fai) file.
    """

    MScompute(f"Gathering length for chromosome {chromosome_name}")

    chrom_lengths = MSpangepopDataHandler.read_fai(fai_file)

    if not chromosome_name.isdigit() or int(chromosome_name) - 1 >= len(chrom_lengths):
        raise MSerror(f"Invalid chromosome number: {chromosome_name}")

    chrom_length = chrom_lengths[int(chromosome_name) - 1]

    if chrom_length <= 0:
        raise MSerror(f"Invalid chromosome length ({chrom_length}) for {chromosome_name}")
    MSsuccess(f"Found chr {chromosome_name} of length: {chrom_length}")

    return chrom_length

def simulate_chromosome_evolution(
    fai_file: str, population_size: int, mutation_rate: float, recombination_rate: float, 
    sample_size: int, output_dir: str, chromosome_name: str, model: str, readable_json, seed
):
    """
    Simulates chromosome evolution and generates a recap file with details.
    """
    pseed = process_seed(seed)
    try:
        start_time = time.time()
        output_dir = os.path.join(output_dir, f"chr_{chromosome_name}")
        os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

        # Fetch chromosome length
        chrom_length = get_chromosome_length(fai_file, chromosome_name)
        MScompute("Starting MSprime simulation...")

        # Create recombination map
        recombination_map = create_recombination_map(chrom_length, recombination_rate)

        # Run ancestry simulation
        ancestry_ts = msprime.sim_ancestry(
            samples=sample_size,
            recombination_rate=recombination_map,
            population_size=population_size,
            random_seed=pseed
        ).simplify()

        # Simulate mutations
        mutated_ts = msprime.sim_mutations(
            ancestry_ts,
            rate=mutation_rate,
            discrete_genome=True,
            model=model,
            random_seed=pseed)

        mutated_ts = mutated_ts.keep_intervals([[0, chrom_length]], simplify=True).trim()

        simulation_time = time.time() - start_time  # Time taken for simulation
        MSsuccess(f"Simulation completed for chromosome {chromosome_name}")

        # Print tree structure (only if there are less than 10 trees)
        total_trees = mutated_ts.num_trees
        if total_trees < 10 and sample_size <= 20:
            print(mutated_ts.draw_text())
                    
            # Save ancestry tree visualization
            ancestry_svg_path = os.path.join(output_dir, f"chr_{chromosome_name}_ancestry.svg")
            with open(ancestry_svg_path, "w") as f:
                f.write(ancestry_ts.draw_svg())

            mutation_svg_path = os.path.join(output_dir, f"chr_{chromosome_name}_mutations.svg")
            with open(mutation_svg_path, "w") as f:
                f.write(mutated_ts.draw_svg())
        else :
            MSsuccess(f"Total number of trees simulated : {total_trees} (see recap file for more data)")

        # Save mutation visualization
        MScompute(f"Saving output for chromosome {chromosome_name}")

        # Record time taken for saving output
        save_start_time = time.time()
        save_output(mutated_ts, chromosome_name, output_dir, readable_json)
        save_time = time.time() - save_start_time  # Time taken to save output

        # Write recap file
        recap_file_path = os.path.join(output_dir, f"chr_{chromosome_name}_simulation_recap.txt")
        with open(recap_file_path, "w") as recap_file:
            recap_file.write("🔹 MSpangepop MSprime Simulation Recap \n")
            recap_file.write("-" * 40 + "\n")
            recap_file.write(f"Seed specified : {seed}\n")
            recap_file.write(f"Chromosome indice (from fai): {chromosome_name}\n")
            recap_file.write(f"Chromosome Length: {chrom_length} bp\n")
            recap_file.write(f"Effective Population Size (Ne): {population_size}\n")
            recap_file.write(f"Mutation Rate: {mutation_rate} per bp\n")
            recap_file.write(f"Recombination Rate: {recombination_rate} per bp\n")
            recap_file.write(f"Sample Size: {sample_size} individuals\n")
            recap_file.write(f"Mutation Model: {model}\n")
            recap_file.write(f"Output Directory: {output_dir}\n\n")
            recap_file.write(f"Total Trees: {mutated_ts.num_trees}\n")
            recap_file.write(f"Total Nodes: {mutated_ts.num_nodes} (can be shared between Trees)\n")
            recap_file.write(f"Total Mutations: {mutated_ts.num_mutations}\n")
            recap_file.write(f"Total Edges: {mutated_ts.num_edges}\n\n")
            recap_file.write(f"Time Taken for Simulation: {simulation_time:.2f} seconds\n")
            recap_file.write(f"Time Taken for Saving Output: {save_time:.2f} seconds\n")
            recap_file.write("-" * 40 + "\n")

        MSsuccess(f"Simulation recap saved to {recap_file_path}")

        # Final timing
        total_time = time.time() - start_time
        MSsuccess(f"Total runtime: {total_time/60:.2f} min.")

    except FileNotFoundError as e:
        raise MSerror(f"MSpangepop -> Missing file: {e}")

    except ValueError as e:
        raise MSerror(f"MSpangepop -> Invalid value encountered: {e}")

    except msprime.LibraryError as e:
        raise MSerror(f"MSpangepop -> msprime simulation error: {e}")

    except Exception as e:
        raise MSerror(f"MSpangepop -> Unexpected error: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulate chromosome evolution using msprime and generate a JSON output."
    )

    # Required arguments with better descriptions
    parser.add_argument('-fai', '--fai', type=str, required=True,
                        help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-p', '--population_size', type=int, required=True,
                        help='Effective population size (Ne). Must be a positive integer.')
    parser.add_argument('-mu', '--mutation_rate', type=float, required=True,
                        help='Mutation rate per base pair (µ). Must be non-negative.')
    parser.add_argument('-r', '--recombination_rate', type=float, required=True,
                        help='Recombination rate per base pair. Must be non-negative.')
    parser.add_argument('-n', '--sample_size', type=int, required=True,
                        help='Sample size (number of individuals to simulate). Must be positive.')
    parser.add_argument('-o', '--output_dir', type=str, required=True,
                        help='Directory to save the output files.')
    parser.add_argument('-c', '--chromosome', type=str, required=True,
                        help='Chromosome number (1-based index). Must be a valid integer string.')
    parser.add_argument('-mo', '--model', type=str, required=True,
                        help='Mutation model to use in msprime simulations.')
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False],
                        help="Save JSON in a human-readable format (True/False, default: False).")
    parser.add_argument('-s', '--seed')
    args = parser.parse_args()

    # Input validation
    if args.population_size <= 0: raise MSerror("Population size must be a positive integer.")
    if args.mutation_rate < 0: raise MSerror("Mutation rate must be non-negative.")
    if args.recombination_rate < 0: raise MSerror("Recombination rate must be non-negative.")
    if args.sample_size <= 0: raise MSerror("Sample size must be a positive integer.")
    if not args.chromosome.isdigit() or int(args.chromosome) <= 0: raise MSerror(f"Invalid chromosome number: {args.chromosome}")
    
    simulate_chromosome_evolution(
        fai_file=args.fai,
        population_size=args.population_size,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        sample_size=args.sample_size,
        output_dir=args.output_dir,
        chromosome_name=args.chromosome,
        model=args.model,
        readable_json=args.readable_json,
        seed=args.seed
    )
