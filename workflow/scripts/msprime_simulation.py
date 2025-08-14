# -*- coding: utf-8 -*-
"""
Author: Sukanya Denni, Lucien Piat
Institution: INRAe
Project: PangenOak

This script generates a JSON file for a specified chromosome using msprime simulations.
It utilizes demographic scenarios from JSON files and outputs the result to a specified directory.
"""

import msprime
import argparse
import os
import json
import time
from collections import defaultdict
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
import matplotlib.pyplot as plt
import demesdraw
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

def load_demographic_model(demographic_file: str, verbose = True):
    """
    Load demographic model from JSON file and create msprime demography object.
    Returns demography object, sample configuration, mutation rate, and recombination rate.
    """
    try:
        with open(demographic_file, 'r') as f:
            demo_data = json.load(f)

        if verbose == True : 
            MScompute(f"Loading demographic model: {demo_data.get('name', 'Unknown')}")

        # Validate required fields
        required_fields = ['populations', 'samples', 'mutation_rate', 'recombination_rate']
        for field in required_fields:
            if field not in demo_data:
                raise MSerror(f"Missing required field '{field}' in demographic JSON file")

        # Create demography object
        demography = msprime.Demography()

        # Add populations
        for pop in demo_data['populations']:
            demography.add_population(
                name=pop['id'],
                initial_size=pop['initial_size'],
                description=pop.get('description', '')
            )

        # Add demographic events
        for event in demo_data.get('demographic_events', []):
            time = event['time']
            event_type = event['type']

            if event_type == 'population_parameters_change':
                demography.add_population_parameters_change(
                    time=time,
                    population=event['population'],
                    initial_size=event.get('size'),
                    growth_rate=event.get('growth_rate', 0)
                )

            elif event_type == 'mass_migration':
                demography.add_mass_migration(
                    time=time,
                    source=event['source'],
                    dest=event['dest'],
                    proportion=event['proportion']
                )

            elif event_type == 'population_split':
                demography.add_population_split(
                    time=time,
                    derived=event['derived'],
                    ancestral=event['ancestral']
                )

            elif event_type == 'add_migration_rate_change':
                demography.add_migration_rate_change(
                    time=time,
                    rate=event['rate'],
                    source=event['source'],
                    dest=event['dest']
                )

        # Set migration rates at time = 0
        if 'migration_matrix' in demo_data and demo_data['migration_matrix']:
            for migration in demo_data['migration_matrix']:
                if migration['time'] == 0:
                    demography.set_migration_rate(
                        source=migration['source'],
                        dest=migration['dest'],
                        rate=migration['rate']
                    )

        # Set up sample configuration
        samples = []
        total_sample_size = 0
        for sample_spec in demo_data['samples']:
            pop_name = sample_spec['population']
            sample_size = sample_spec['sample_size']
            samples.append(msprime.SampleSet(sample_size, population=pop_name))
            total_sample_size += sample_size

        mutation_rate = demo_data['mutation_rate']
        recombination_rate = demo_data['recombination_rate']

        if mutation_rate < 0:
            raise MSerror("Mutation rate must be non-negative.")
        if recombination_rate < 0:
            raise MSerror("Recombination rate must be non-negative.")

        demography.sort_events()

        return demography, samples, total_sample_size, mutation_rate, recombination_rate, demo_data

    except FileNotFoundError:
        raise MSerror(f"Demographic file not found: {demographic_file}")
    except json.JSONDecodeError as e:
        raise MSerror(f"Invalid JSON in demographic file: {e}")
    except KeyError as e:
        raise MSerror(f"Missing key in demographic file: {e}")
    except Exception as e:
        raise MSerror(f"Error loading demographic model: {e}")

def plot_demographic(demographic_file: str, output_file: str, log_time: bool = True):
    """
    Plot demographic model using demes visualization.
    
    Args:
        demographic_file: Path to demographic JSON file
        output_file: Path to save the output plot
        log_time: Whether to use log scale for time axis
    """
        
    # Load the demographic model
    demography, _, _, _, _, demo_data = load_demographic_model(demographic_file, verbose = False)
    
    # Create figure and plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Convert to demes and draw
    graph = demography.to_demes()
    demesdraw.tubes(graph, ax=ax, log_time=log_time)
    ax.set_title(f"{demo_data.get('name', 'Demographic Model')}\n{demo_data.get('description', '')}")
    
    # Save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    

def save_output(mutated_ts, chromosome_name: str, json_file: str, readable_json: bool):
    """
    Optimized version with pre-indexed mutations
    O(M + T × N) M is no more a factor
    """
    
    json_filename = json_file
    
    # First, collect all tree intervals and their indices
    tree_intervals = []
    for tree in mutated_ts.trees():
        tree_intervals.append((tree.interval.left, tree.interval.right, tree.index))
    
    # Pre-index all mutations by tree - O(M) 
    tree_mutations = defaultdict(list)
    for mutation in mutated_ts.mutations():
        position = mutated_ts.site(mutation.site).position
        
        # Find which tree this mutation belongs to
        for left, right, tree_idx in tree_intervals:
            if left <= position < right:
                tree_mutations[tree_idx].append({"node": mutation.node})
                break  # Each mutation belongs to only one tree
    
    try:
        with open(json_filename, 'w') as f:
            f.write('[\n')
            first_tree = True
            
            for tree in mutated_ts.trees():
                interval = [int(tree.interval[0]), int(tree.interval[1])]
                
                # Build nodes and edges as before
                nodes = [{"id": node, "time": mutated_ts.node(node).time} for node in tree.nodes()]
                edges = []
                for parent in tree.nodes():
                    children = list(tree.children(parent))
                    edges.extend({"parent": parent, "child": child} for child in children)
                
                # Get mutations for this specific tree - O(1) lookup!
                mutations = tree_mutations.get(tree.index, [])
                
                tree_data = {
                    "tree_index": tree.index,
                    "interval": interval,
                    "nodes": nodes,
                    "edges": edges,
                    "mutations": mutations,
                }
                
                if not first_tree:
                    f.write(",\n")
                first_tree = False
                
                json.dump(tree_data, f, indent=get_indent(readable_json))
            
            f.write('\n]')
    
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

    return chrom_length

def save_tree_sequences(mutated_ts, ancestry_ts, full_ancestry_ts, mutated_ts_file, ancestry_ts_file, full_ancestry_ts_file):
    """Save both tree sequences to files"""
    # Save as .trees files (tskit native format - recommended)
    mutated_ts.dump(mutated_ts_file)
    ancestry_ts.dump(ancestry_ts_file)
    full_ancestry_ts.dump(full_ancestry_ts_file)

def simulate_chromosome_evolution(
    fai_file: str, json_file: str, chromosome_name: str, model: str, readable_json, 
    seed, mutated_ts_file, ancestry_ts_file, full_ancestry_ts_file, recap, demographic_file: str
):
    """
    Simulates chromosome evolution using demographic scenarios from JSON files.
    """
    pseed = process_seed(seed)
    try:
        start_time = time.time()

        # Fetch chromosome length
        chrom_length = get_chromosome_length(fai_file, chromosome_name)
        MScompute(f"Chr {chromosome_name}, starting MSprime simulation...")

        # Load demographic model
        demography, samples, sample_size, mutation_rate, recombination_rate, demo_data = load_demographic_model(demographic_file)
        plot_demographic(demographic_file, recap.replace('.txt', '_plot.png'))
        MScompute(f"Using mutation rate: {mutation_rate}")
        MScompute(f"Using recombination rate: {recombination_rate}")
        MScompute(f"Total sample size: {sample_size}")

        # Create recombination map
        recombination_map = create_recombination_map(chrom_length, recombination_rate)

        # Run ancestry simulation
        full_ancestry_ts = msprime.sim_ancestry(
            samples=samples,
            demography=demography,
            recombination_rate=recombination_map,
            random_seed=pseed, 
            record_full_arg=True
        )
        ancestry_ts = full_ancestry_ts.simplify()
        # Simulate mutations
        mutated_ts = msprime.sim_mutations(
            ancestry_ts,
            rate=mutation_rate,
            discrete_genome=True,
            model=model,
            random_seed=pseed)

        mutated_ts = mutated_ts.simplify()

        simulation_time = time.time() - start_time  # Time taken for simulation

        # Save mutation visualization
        MScompute(f"Saving output for chromosome {chromosome_name}")

        # Record time taken for saving output
        save_start_time = time.time()
        save_output(mutated_ts, chromosome_name, json_file, readable_json)
        save_time = time.time() - save_start_time 

        save_tree_sequences(mutated_ts, ancestry_ts, full_ancestry_ts, mutated_ts_file, ancestry_ts_file, full_ancestry_ts_file)

        # Write recap file
        recap_file_path = recap
        with open(recap_file_path, "w") as recap_file:
            recap_file.write("🔹 MSpangepop MSprime Simulation Recap \n")
            recap_file.write("-" * 40 + "\n")
            recap_file.write(f"Demographic Model: {demo_data.get('name', 'Unknown')}\n")
            recap_file.write(f"Demographic File: {demographic_file}\n")
            recap_file.write(f"Seed specified: {seed}\n")
            recap_file.write(f"Chromosome indice (from fai): {chromosome_name}\n")
            recap_file.write(f"Chromosome Length: {chrom_length} bp\n")
            recap_file.write(f"Mutation Rate: {mutation_rate} per bp\n")
            recap_file.write(f"Recombination Rate: {recombination_rate} per bp\n")
            recap_file.write(f"Sample Size: {sample_size} individuals\n")
            recap_file.write(f"Mutation Model: {model}\n")
            recap_file.write(f"Total Trees: {mutated_ts.num_trees}\n")
            recap_file.write(f"Total Nodes: {mutated_ts.num_nodes} (can be shared between Trees)\n")
            recap_file.write(f"Total Mutations: {mutated_ts.num_mutations}\n")
            recap_file.write(f"Total Edges: {mutated_ts.num_edges}\n\n")
            
            # Population details
            recap_file.write("Population Configuration:\n")
            for pop in demo_data['populations']:
                recap_file.write(f"  - {pop['id']}: {pop.get('description', 'No description')}\n")
                recap_file.write(f"    Initial size: {pop['initial_size']}\n")
            
            recap_file.write("\nSample Distribution:\n")
            for sample in demo_data['samples']:
                recap_file.write(f"  - {sample['population']}: {sample['sample_size']} samples\n")
            
            recap_file.write(f"\nTime Taken for Simulation: {simulation_time:.2f} seconds\n")
            recap_file.write(f"Time Taken for Saving Output: {save_time:.2f} seconds\n")
            recap_file.write("-" * 40 + "\n")
        
        total_time = time.time() - start_time
        MSsuccess(f"Chr {chromosome_name} msprime recap created, total runtime: {total_time/60:.2f} min.")

    except FileNotFoundError as e:
        raise MSerror(f"MSpangepop -> Missing file: {e}")

    except ValueError as e:
        raise MSerror(f"MSpangepop -> Invalid value encountered: {e}")

    except Exception as e:
        raise MSerror(f"MSpangepop -> Unexpected error: {e}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Simulate chromosome evolution using msprime with demographic scenarios from JSON files."
    )

    # Required arguments
    parser.add_argument('-fai', '--fai', type=str, required=True,
                        help='Path to the FAI index file for the reference FASTA.')
    parser.add_argument('-d', '--demographic', type=str, required=True,
                        help='Path to JSON file containing demographic scenario (required).')
    parser.add_argument('--json', type=str, required=True,
                        help='Path to save the output JSON file.')
    parser.add_argument('-c', '--chromosome', type=str, required=True,
                        help='Chromosome number (1-based index). Must be a valid integer string.')
    parser.add_argument('-mo', '--model', type=str, required=True,
                        help='Mutation model to use in msprime simulations (e.g., JC69, HKY, GTR).')
    
    # Optional arguments
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False], default=False,
                        help="Save JSON in a human-readable format (True/False, default: False).")
    parser.add_argument('-s', '--seed', type=str,
                        help='Random seed for reproducibility.')
    parser.add_argument('--mutated_ts', type=str, required=True,
                        help='Path to save the mutated tree sequence file.')
    parser.add_argument('--ancestry_ts', type=str, required=True,
                        help='Path to save the ancestry tree sequence file.')
    parser.add_argument('--full_ancestry_ts', type=str, required=True,
                        help='Path to save the ancestry tree sequence file.')
    parser.add_argument('--recap', type=str, required=True,
                        help='Path to save the simulation recap file.')
    
    args = parser.parse_args()

    # Input validation
    if not os.path.exists(args.demographic):
        raise MSerror(f"Demographic file not found: {args.demographic}")
    
    if not args.chromosome.isdigit() or int(args.chromosome) <= 0:
        raise MSerror(f"Invalid chromosome number: {args.chromosome}")
    
    simulate_chromosome_evolution(
        fai_file=args.fai,
        json_file=args.json,
        chromosome_name=args.chromosome,
        model=args.model,
        readable_json=args.readable_json,
        seed=args.seed,
        mutated_ts_file=args.mutated_ts,
        ancestry_ts_file=args.ancestry_ts,
        full_ancestry_ts_file=args.full_ancestry_ts,
        recap=args.recap,
        demographic_file=args.demographic
    )