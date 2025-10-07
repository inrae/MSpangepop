"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Augment JSON file with variant type and size.
"""

import argparse
import random
import os
import time
import json
import concurrent.futures
from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute, MSwarning, process_seed

class Selector:
    """Static methods for selecting mutation attributes."""
    
    @staticmethod
    def position(interval) -> int:
        """Returns a random position within the given interval."""
        return random.randint(interval[0] , interval[1]-2)

    @staticmethod
    def type(variant_probabilities) -> str:
        """Selects a variant type based on predefined probabilities."""
        try:
            # Normalize probabilities to ensure they sum to 1
            total = sum(variant_probabilities.values())
            if total == 0:
                raise MSerror("All variant probabilities are zero")
            
            normalized_probs = {k: v/total for k, v in variant_probabilities.items()}
            return random.choices(
                list(normalized_probs.keys()), 
                weights=list(normalized_probs.values())
            )[0]
        except Exception as e:
            raise MSerror(f"Error in selecting variant type: {e}")

    @staticmethod 
    def length(length_df, max_length, minimal_sv_length):
        try:
            rand_val = random.random()
            row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]
            lower_bound, upper_bound = map(float, row['size_interval'].strip('[]').split(','))
            length = random.randint(int(lower_bound), int(upper_bound))
            if max_length:
                return max(min(minimal_sv_length, max_length), min(length, max_length))
            else:
                return max(minimal_sv_length, length)
        except Exception as e:
            raise MSerror(f"Error in selecting variant length: {e}")

def validate_sv_distribution(sv_dist):
    """Validate that SV distribution is properly formatted and sums to 100 (or can be normalized)."""
    required_types = {'SNP', 'DEL', 'INS', 'INV', 'DUP'}
    
    # Check all required types are present
    missing_types = required_types - set(sv_dist.keys())
    if missing_types:
        # Add missing types with 0 probability
        for sv_type in missing_types:
            sv_dist[sv_type] = 0
    
    # Check that at least one type has non-zero probability
    if sum(sv_dist.values()) == 0:
        raise MSerror("All SV type probabilities are zero")
    
    # Warn if not summing to 100 (will be normalized anyway)
    total = sum(sv_dist.values())
    if abs(total - 100) > 0.01:
        MSwarning(f"SV distribution sums to {total}, will be normalized to 100%")
    
    return sv_dist

def augment_type(tree, variant_probabilities):
    """Augments mutations in a tree with variant type"""
    try:
        for node in tree.get("nodes", []):
            for mutation in node.get("mutations", []):
                # Assign variant type
                mutation["type"] = Selector.type(variant_probabilities)
        return tree
    except Exception as e:
        MSerror(f"Error augmenting mutations: {e}")
        return None

def augment_length_and_position(tree, length_files, minimal_sv_length):
    """Augments mutations in a tree with length and position based on variant types.
    
    Simple boundary protection: just ensure we don't place mutations at obvious
    boundary positions. The actual validation happens in apply_mutations_to_graphs.
    """
    fails = 0
    boundary_fails = 0
    
    try:
        for node in tree.get("nodes", []):  
            for mutation in node.get("mutations", []):
                
                interval_start = node["interval"][0]
                interval_end = node["interval"][1]
                interval_size = interval_end - interval_start
                
                # Check if interval is too small for any mutation
                if interval_size < 3:
                    fails += 1
                    mutation["type"] = None
                    mutation["length"] = None
                    mutation["start"] = None
                    continue

                variant_type = mutation["type"]
                
                # For SNPs: Avoid first and last position
                if variant_type == "SNP":
                    if interval_size <= 2:  # No valid interior positions
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    # Select position excluding first and last base
                    start_pos = random.randint(interval_start + 1, interval_end - 2)
                    length = None
                    
                # For INS: Avoid first position only
                elif variant_type == "INS":
                    # Can insert anywhere except at position 0
                    start_pos = random.randint(interval_start + 1, interval_end - 1)
                    
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, None, minimal_sv_length)
                    
                    # Update affected node intervals for tracking
                    affected_nodes = node.get("affected_nodes", [])
                    for affected_node_id in affected_nodes:
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] += length
                    
                # For DEL: Ensure we don't delete first/last and leave >=3 nodes
                elif variant_type == "DEL":
                    # Need room to delete something without affecting boundaries
                    if interval_size <= 3:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    # Start after first position, leave room for at least one deletion
                    max_start = interval_end - 2  # Need at least 1 position to delete
                    if interval_start + 1 > max_start:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    start_pos = random.randint(interval_start + 1, max_start)
                    
                    # Maximum deletion: don't delete the last position
                    max_del_length = interval_end - 1 - start_pos
                    
                    # Also ensure we keep interval >= 3
                    max_del_length = min(max_del_length, interval_size - 3)
                    
                    if max_del_length < 1:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, max_del_length, 1)
                    
                    # Update affected node intervals for tracking
                    affected_nodes = node.get("affected_nodes", [])
                    for affected_node_id in affected_nodes:
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] -= length
                
                # For DUP and INV: Don't affect first or last position
                elif variant_type in ["DUP", "INV"]:
                    # Need at least 2 positions to duplicate/invert (can't be first or last)
                    if interval_size <= 3:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    # Start after first, need room for at least minimal_sv_length
                    max_start = interval_end - 1 - minimal_sv_length
                    
                    if interval_start + 1 > max_start:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    start_pos = random.randint(interval_start + 1, max_start)
                    
                    # Maximum length: don't include the last position
                    max_length = interval_end - 1 - start_pos
                    
                    if max_length < minimal_sv_length:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, max_length, minimal_sv_length)
                    
                    # For DUP, update affected node intervals for tracking
                    if variant_type == "DUP":
                        affected_nodes = node.get("affected_nodes", [])
                        for affected_node_id in affected_nodes:
                            affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                            if affected_node:
                                affected_node["interval"][1] += length
                
                else:
                    # Unknown variant type
                    mutation["type"] = None
                    mutation["length"] = None
                    mutation["start"] = None
                    continue
                
                # Assign the calculated values
                mutation["length"] = length
                mutation["start"] = start_pos

        return tree, fails, boundary_fails
        
    except Exception as e:
        raise MSerror(f"Error augmenting mutations with length and position: {e}")

def process_single_tree(tree, variant_probabilities, minimal_sv_length, seed_offset):
    """
    Processes a single tree. Loads length files inside worker to avoid pickling.
    """
    random.seed(seed_offset)
    
    # Load length files inside worker (avoid pickling large dataframes)
    length_files = {}
    for var_type, file_path in {
        'DEL': 'simulation_data/size_distribDEL.tsv',
        'INS': 'simulation_data/size_distribINS.tsv',
        'INV': 'simulation_data/size_distribINV.tsv',
        'DUP': 'simulation_data/size_distribDUP.tsv'
    }.items():
        if os.path.exists(file_path):
            length_files[var_type] = MSpangepopDataHandler.read_variant_length_file(file_path)
    
    try:
        augmented_tree = augment_type(tree, variant_probabilities)
        augmented_tree, fails, boundary_fails = augment_length_and_position(
            augmented_tree, length_files, minimal_sv_length
        )
        return augmented_tree, fails, boundary_fails
    except MSerror as e:
        raise MSerror(f"Error processing the tree: {e}")

def process_batch(batch, variant_probabilities, num_threads, minimal_sv_length, seed_offset_start):
    """Process a batch of trees in parallel."""
    if seed_offset_start is None:
        # Non-deterministic: generate random seeds
        seed_list = [random.randint(0, 2**32 - 1) for _ in range(len(batch))]
    else:
        # Deterministic: use sequential seeds
        seed_list = [seed_offset_start + i for i in range(len(batch))]
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(
            process_single_tree, 
            batch,
            [variant_probabilities] * len(batch), 
            [minimal_sv_length] * len(batch),
            seed_list
        ))
    return results
    
def stream_json_trees(json_file):
    """Generator that yields trees one at a time from JSON array."""
    with open(json_file, 'r') as f:
        # Skip opening bracket
        content = f.read().strip()
        if content.startswith('['):
            content = content[1:]
        if content.endswith(']'):
            content = content[:-1]
        
        # Split by top-level objects (simple approach)
        # For production, use ijson for proper streaming
        trees = json.loads('[' + content + ']')
        for tree in trees:
            yield tree

def main(json_file, output_json_file, sv_distribution, chromosome, num_threads, 
         minimal_sv_length, readable_json, seed, batch_size=1000):
    """Processes trees in batches to handle millions of trees efficiently."""
    try:
        start_time = time.time()
        MScompute(f"Generating variants for chromosome {chromosome} using {num_threads} threads")
        MScompute(f"Processing in batches of {batch_size} trees")
        
        # Parse and validate SV distribution
        if isinstance(sv_distribution, str):
            try:
                variant_probabilities = json.loads(sv_distribution)
            except json.JSONDecodeError as e:
                raise MSerror(f"Invalid JSON format for SV distribution: {e}")
        else:
            variant_probabilities = sv_distribution
            
        variant_probabilities = validate_sv_distribution(variant_probabilities)
        
        sv_summary = ", ".join([f"{k}:{v}%" for k, v in variant_probabilities.items()])
        MScompute(f"Using SV distribution: {sv_summary}")

        # Set deterministic seed (can be None for non-deterministic)
        seed = process_seed(seed)
        if seed is not None:
            random.seed(seed)

        # Initialize counters
        mutation_counts = {sv_type: 0 for sv_type in variant_probabilities.keys()}
        total_mutations = 0
        total_trees = 0
        total_fails = 0
        total_boundary_fails = 0
        
        # Open output file and write opening bracket
        with open(output_json_file, 'w') as outfile:
            outfile.write('[\n' if readable_json else '[')
            
            first_batch = True
            batch = []
            tree_index = 0
            
            # Stream and process trees in batches
            for tree in stream_json_trees(json_file):
                batch.append(tree)
                tree_index += 1
                
                if len(batch) >= batch_size:
                    # Process batch
                    # Handle None seed case
                    if seed is None:
                        seed_offset = None
                    else:
                        seed_offset = seed + total_trees
                    
                    results = process_batch(batch, variant_probabilities, num_threads, 
                                          minimal_sv_length, seed_offset)
                    
                    # Write results and update counts
                    for augmented_tree, fails, boundary_fails in results:
                        if not first_batch:
                            outfile.write(',\n' if readable_json else ',')
                        else:
                            first_batch = False
                        
                        json.dump(augmented_tree, outfile, indent=2 if readable_json else None)
                        
                        # Update counters
                        total_fails += fails
                        total_boundary_fails += boundary_fails
                        total_trees += 1
                        
                        for node in augmented_tree.get("nodes", []):
                            for mutation in node.get("mutations", []):
                                if mutation.get("type"):
                                    mutation_counts[mutation["type"]] += 1
                                    total_mutations += 1
                    
                    MScompute(f"Processed {total_trees} trees...")
                    batch = []
            
            # Process remaining trees in last batch
            if batch:
                # Handle None seed case
                if seed is None:
                    seed_offset = None
                else:
                    seed_offset = seed + total_trees
                
                results = process_batch(batch, variant_probabilities, num_threads, 
                                      minimal_sv_length, seed_offset)
                
                for augmented_tree, fails, boundary_fails in results:
                    if not first_batch:
                        outfile.write(',\n' if readable_json else ',')
                    else:
                        first_batch = False
                    
                    json.dump(augmented_tree, outfile, indent=2 if readable_json else None)
                    
                    total_fails += fails
                    total_boundary_fails += boundary_fails
                    total_trees += 1
                    
                    for node in augmented_tree.get("nodes", []):
                        for mutation in node.get("mutations", []):
                            if mutation.get("type"):
                                mutation_counts[mutation["type"]] += 1
                                total_mutations += 1
            
            # Close JSON array
            outfile.write('\n]' if readable_json else ']')

        # Print summary
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        if total_fails > 0:
            MSwarning(f"{total_fails} mutations skipped due to intervals < 3 bases")
        if total_boundary_fails > 0:
            MSwarning(f"{total_boundary_fails} mutations rejected due to boundary constraints")
        
        MScompute(f"Total trees processed: {total_trees}")
        MScompute(f"Total mutations: {total_mutations}")
        MScompute(f"Mutation breakdown: {', '.join([f'{k}:{v}' for k, v in mutation_counts.items()])}")
        MSsuccess(f"Processing time: {elapsed_time:.2f} seconds ({total_trees/elapsed_time:.1f} trees/sec)")
        

    except Exception as e:
        raise MSerror(f"Critical error processing (Chromosome {chromosome}): {e}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--sv_distribution", required=True, 
                        help="SV type distribution as JSON string")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--minimal_sv_length", type=int, default=1, help="Minimal size for variants generated")
    parser.add_argument("--batch_size", type=int, default=1000, help="Number of trees to process per batch")
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False], default=False,
                        help="Save JSON in a human-readable format (True/False, default: False).")
    parser.add_argument("--seed", help="Random seed for reproducibility")
    args = parser.parse_args()

    if args.minimal_sv_length < 1:
        raise MSerror('minimal_sv_length cant be less than 1')

    main(args.json, args.output, args.sv_distribution, args.chromosome, args.threads, 
         args.minimal_sv_length, args.readable_json, args.seed, args.batch_size)