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

        if fails != 0:
            MSwarning(f"{fails} mutations skipped due to intervals < 3 bases")
        if boundary_fails != 0:
            MSwarning(f"{boundary_fails} mutations rejected due to boundary constraints")
        
        return tree
        
    except Exception as e:
        raise MSerror(f"Error augmenting mutations with length and position: {e}")

def process_single_tree(tree, length_files, variant_probabilities, minimal_sv_length, seed_offset):
    """
    Processes a single tree, augmenting mutation types, positions, and lengths.
    A deterministic seed is set here so each tree gets reproducible random results
    regardless of which process executes it.
    """
    random.seed(seed_offset)  # <-- NEW: Set deterministic seed per tree for reproducibility
    try:
        # First augment with mutation types
        augmented_tree = augment_type(tree, variant_probabilities)
        
        # Now augment with lengths and positions
        augmented_tree = augment_length_and_position(augmented_tree, length_files, minimal_sv_length)
        
        return augmented_tree
    except MSerror as e:
        raise MSerror(f"Error processing the tree: {e}")

def process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_sv_length, base_seed):
    """
    Parallel processing of trees, passing a unique deterministic seed to each worker.
    If base_seed is None, seeds are generated randomly for non-deterministic runs.
    """
    if base_seed is None:
        # Generate a fresh random seed for each tree (non-deterministic)
        seed_list = [random.randint(0, 2**32 - 1) for _ in range(len(tree_list))]
    else:
        # Deterministic: base_seed + index
        seed_list = [base_seed + i for i in range(len(tree_list))]

    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(
            process_single_tree, 
            tree_list, 
            [length_files] * len(tree_list), 
            [variant_probabilities] * len(tree_list), 
            [minimal_sv_length] * len(tree_list),
            seed_list
        ))
    return results


def main(json_file, output_json_file, sv_distribution, chromosome, num_threads, minimal_sv_length, readable_json, seed):
    """Processes a JSON list of trees, augmenting mutations with variant type and size."""
    try:
        start_time = time.time()
        MScompute(f"Generating variants for chromosome {chromosome} using {num_threads} threads")
        
        # Parse and validate SV distribution
        if isinstance(sv_distribution, str):
            try:
                variant_probabilities = json.loads(sv_distribution)
            except json.JSONDecodeError as e:
                raise MSerror(f"Invalid JSON format for SV distribution: {e}")
        else:
            variant_probabilities = sv_distribution
            
        variant_probabilities = validate_sv_distribution(variant_probabilities)
        
        # Log the distribution being used
        sv_summary = ", ".join([f"{k}:{v}%" for k, v in variant_probabilities.items()])
        MScompute(f"Using SV distribution: {sv_summary}")

        # Set deterministic seed for the main process
        seed = process_seed(seed)
        random.seed(seed)

        # Read tree data (only once)
        tree_list = MSpangepopDataHandler.read_json(json_file)
        if not isinstance(tree_list, list) or not tree_list:
            raise MSerror('Tree JSON file is empty or improperly formatted.')

        # Read length distribution files (only once)
        length_files = {}
        for var_type, file_path in {
            'DEL': 'simulation_data/size_distribDEL.tsv',
            'INS': 'simulation_data/size_distribINS.tsv',
            'INV': 'simulation_data/size_distribINV.tsv',
            'DUP': 'simulation_data/size_distribDUP.tsv'
        }.items():
            if not os.path.exists(file_path):
                MSwarning(f"Length distribution file missing for {var_type}")
            else:
                length_files[var_type] = MSpangepopDataHandler.read_variant_length_file(file_path)

        # Process and augment mutations in parallel with reproducible seeds
        tree_list = process_trees_in_parallel(tree_list, length_files, variant_probabilities, num_threads, minimal_sv_length, seed)

        # Count total mutations by type
        mutation_counts = {sv_type: 0 for sv_type in variant_probabilities.keys()}
        total_mutations = 0
        
        for tree in tree_list:
            for node in tree.get("nodes", []):
                for mutation in node.get("mutations", []):
                    if mutation.get("type"):
                        mutation_counts[mutation["type"]] += 1
                        total_mutations += 1
        
        # Save output
        MSpangepopDataHandler.save_json(tree_list, output_json_file, readable_json=readable_json)

        # Print summary
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        
        MScompute(f"Total mutations: {total_mutations}")
        MScompute(f"Mutation breakdown: {', '.join([f'{k}:{v}' for k, v in mutation_counts.items()])}")
        MScompute(f"Processing time: {elapsed_time:.2f} seconds")
        

    except Exception as e:
        raise MSerror(f"Critical error processing (Chromosome {chromosome}): {e}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Augment JSON file with variant type and size.")
    parser.add_argument("--json", required=True, help="Path to the JSON file containing tree and mutation data.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file where augmented data will be saved.")
    parser.add_argument("--sv_distribution", required=True, 
                        help="SV type distribution as JSON string, e.g., '{\"SNP\": 50, \"DEL\": 20, \"INS\": 20, \"INV\": 10, \"DUP\": 0}'")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--minimal_sv_length", type=int, default=1, help="Minimal size for variants generated")
    parser.add_argument("--readable_json", type=lambda x: x.lower() == 'true',
                        choices=[True, False], default=False,
                        help="Save JSON in a human-readable format (True/False, default: False).")
    parser.add_argument("--seed", help="Random seed for reproducibility")
    args = parser.parse_args()

    if args.minimal_sv_length < 1:
        raise MSerror('minimal_sv_length cant be less than 1')

    main(args.json, args.output, args.sv_distribution, args.chromosome, args.threads, args.minimal_sv_length, args.readable_json, args.seed)
