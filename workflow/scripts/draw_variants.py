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
        return random.randint(interval[0], interval[1]-2)

    @staticmethod
    def type(variant_probabilities) -> str:
        """Selects a variant type based on predefined probabilities."""
        try:
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
    def type_with_max(variant_probabilities, current_counts, max_sv):
        """
        Selects a variant type respecting max_sv limits.
        Returns None if all types have reached their max.
        """
        available_types = {
            k: v for k, v in variant_probabilities.items()
            if max_sv.get(k, float('inf')) > current_counts.get(k, 0)
        }
        
        if not available_types:
            return None
        
        total = sum(available_types.values())
        if total == 0:
            return random.choice(list(available_types.keys()))
        
        normalized_probs = {k: v/total for k, v in available_types.items()}
        return random.choices(
            list(normalized_probs.keys()), 
            weights=list(normalized_probs.values())
        )[0]

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


def load_simulation_params(demographic_file):
    """Load simulation parameters from demographic JSON file."""
    with open(demographic_file, 'r') as f:
        demo_data = json.load(f)
    
    sim_params = demo_data.get('simulation_params', {})
    
    return {
        'sv_distribution': sim_params.get('sv_distribution', {}),
        'max_sv': sim_params.get('max_sv'),
        'sv_length_files': sim_params.get('sv_length_files'),
        'minimal_sv_length': sim_params.get('minimal_sv_length', 1),
        'readable_json': sim_params.get('readable_json', False),
        'seed': sim_params.get('seed')
    }


def validate_sv_distribution(sv_dist):
    """Validate that SV distribution is properly formatted."""
    required_types = {'SNP', 'DEL', 'INS', 'INV', 'DUP'}
    
    missing_types = required_types - set(sv_dist.keys())
    for sv_type in missing_types:
        sv_dist[sv_type] = 0
    
    if sum(sv_dist.values()) == 0:
        raise MSerror("All SV type probabilities are zero")
    
    total = sum(sv_dist.values())
    if abs(total - 100) > 0.01:
        MSwarning(f"SV distribution sums to {total}, will be normalized to 100%")
    
    return sv_dist


def validate_max_sv(max_sv):
    """Validate max_sv parameter."""
    if max_sv is None:
        return None
    
    required_types = {'SNP', 'DEL', 'INS', 'INV', 'DUP'}
    
    for sv_type in required_types:
        if sv_type not in max_sv:
            max_sv[sv_type] = float('inf')
        elif max_sv[sv_type] is None:
            max_sv[sv_type] = float('inf')
    
    for sv_type in max_sv:
        if max_sv[sv_type] != float('inf'):
            max_sv[sv_type] = int(max_sv[sv_type])
    
    limited = {k: v for k, v in max_sv.items() if v != float('inf')}
    if limited:
        MScompute(f"Max SV limits: {', '.join([f'{k}:{v}' for k, v in limited.items()])}")
        MScompute(f"Total max mutations: {sum(limited.values())}")
    
    return max_sv


def load_length_files_worker(sv_length_files):
    """Load length files inside worker process."""
    default_files = {
        'DEL': 'simulation_data/size_distribDEL.tsv',
        'INS': 'simulation_data/size_distribINS.tsv',
        'INV': 'simulation_data/size_distribINV.tsv',
        'DUP': 'simulation_data/size_distribDUP.tsv'
    }
    
    file_paths = default_files.copy()
    if sv_length_files:
        file_paths.update(sv_length_files)
    
    length_files = {}
    for var_type, file_path in file_paths.items():
        if file_path and os.path.exists(file_path):
            length_files[var_type] = MSpangepopDataHandler.read_variant_length_file(file_path)
    
    return length_files


def stream_json_trees(json_file):
    """Generator that yields trees one at a time from JSON array."""
    with open(json_file, 'r') as f:
        content = f.read().strip()
        if content.startswith('['):
            content = content[1:]
        if content.endswith(']'):
            content = content[:-1]
        
        trees = json.loads('[' + content + ']')
        for tree in trees:
            yield tree


def collect_all_mutations(json_file):
    """
    First pass: Collect all mutation addresses from all trees.
    Returns list of tuples: (tree_idx, node_idx, mutation_idx, interval_size)
    """
    MScompute("First pass: Collecting all mutation positions...")
    
    mutation_addresses = []
    tree_idx = 0
    
    for tree in stream_json_trees(json_file):
        for node_idx, node in enumerate(tree.get("nodes", [])):
            interval_start = node["interval"][0]
            interval_end = node["interval"][1]
            interval_size = interval_end - interval_start
            
            for mut_idx, _ in enumerate(node.get("mutations", [])):
                mutation_addresses.append((tree_idx, node_idx, mut_idx, interval_size))
        
        tree_idx += 1
        if tree_idx % 10000 == 0:
            MScompute(f"  Scanned {tree_idx} trees, found {len(mutation_addresses)} mutations...")
    
    MScompute(f"Total mutations found: {len(mutation_addresses)}")
    return mutation_addresses, tree_idx


def assign_global_types(mutation_addresses, variant_probabilities, max_sv, seed):
    """
    Assign variant types globally, respecting max_sv limits.
    Returns dict mapping (tree_idx, node_idx, mut_idx) -> variant_type (or None if deleted)
    """
    if seed is not None:
        random.seed(seed)
    
    # Shuffle to avoid bias toward early trees
    shuffled_addresses = mutation_addresses.copy()
    random.shuffle(shuffled_addresses)
    
    assignments = {}
    current_counts = {sv_type: 0 for sv_type in variant_probabilities.keys()}
    
    total_max = sum(v for v in max_sv.values() if v != float('inf'))
    has_limits = total_max > 0
    
    selected_count = 0
    
    for tree_idx, node_idx, mut_idx, interval_size in shuffled_addresses:
        if interval_size < 3:
            assignments[(tree_idx, node_idx, mut_idx)] = None
            continue
        
        if has_limits and selected_count >= total_max:
            assignments[(tree_idx, node_idx, mut_idx)] = None
            continue
        
        variant_type = Selector.type_with_max(variant_probabilities, current_counts, max_sv)
        
        if variant_type is None:
            assignments[(tree_idx, node_idx, mut_idx)] = None
        else:
            assignments[(tree_idx, node_idx, mut_idx)] = variant_type
            current_counts[variant_type] += 1
            selected_count += 1
    
    MScompute(f"Assigned {selected_count} mutations: {current_counts}")
    deleted = len(mutation_addresses) - selected_count
    if deleted > 0:
        MScompute(f"Mutations deleted (exceeded limits): {deleted}")
    
    return assignments, current_counts


def augment_type(tree, variant_probabilities):
    """Augments mutations in a tree with variant type (distribution mode)."""
    try:
        for node in tree.get("nodes", []):
            for mutation in node.get("mutations", []):
                mutation["type"] = Selector.type(variant_probabilities)
        return tree
    except Exception as e:
        MSerror(f"Error augmenting mutations: {e}")
        return None


def augment_type_from_assignments(tree, tree_idx, assignments):
    """Augments mutations in a tree using pre-computed assignments (max_sv mode)."""
    try:
        for node_idx, node in enumerate(tree.get("nodes", [])):
            mutations_to_keep = []
            for mut_idx, mutation in enumerate(node.get("mutations", [])):
                assigned_type = assignments.get((tree_idx, node_idx, mut_idx))
                if assigned_type is not None:
                    mutation["type"] = assigned_type
                    mutations_to_keep.append(mutation)
            node["mutations"] = mutations_to_keep
        return tree
    except Exception as e:
        MSerror(f"Error augmenting mutations from assignments: {e}")
        return None


def augment_length_and_position(tree, length_files, minimal_sv_length):
    """Augments mutations in a tree with length and position based on variant types."""
    fails = 0
    boundary_fails = 0
    
    try:
        for node in tree.get("nodes", []):  
            for mutation in node.get("mutations", []):
                
                interval_start = node["interval"][0]
                interval_end = node["interval"][1]
                interval_size = interval_end - interval_start
                
                if interval_size < 3:
                    fails += 1
                    mutation["type"] = None
                    mutation["length"] = None
                    mutation["start"] = None
                    continue

                variant_type = mutation["type"]
                
                if variant_type == "SNP":
                    if interval_size <= 2:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    start_pos = random.randint(interval_start + 1, interval_end - 2)
                    length = None
                    
                elif variant_type == "INS":
                    start_pos = random.randint(interval_start + 1, interval_end - 1)
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, None, minimal_sv_length) if length_df is not None else minimal_sv_length
                    
                    for affected_node_id in node.get("affected_nodes", []):
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] += length
                    
                elif variant_type == "DEL":
                    if interval_size <= 3:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    max_start = interval_end - 2
                    if interval_start + 1 > max_start:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    start_pos = random.randint(interval_start + 1, max_start)
                    max_del_length = min(interval_end - 1 - start_pos, interval_size - 3)
                    
                    if max_del_length < 1:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, max_del_length, 1) if length_df is not None else 1
                    
                    for affected_node_id in node.get("affected_nodes", []):
                        affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                        if affected_node:
                            affected_node["interval"][1] -= length
                
                elif variant_type in ["DUP", "INV"]:
                    if interval_size <= 3:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    max_start = interval_end - 1 - minimal_sv_length
                    if interval_start + 1 > max_start:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    start_pos = random.randint(interval_start + 1, max_start)
                    max_length = interval_end - 1 - start_pos
                    
                    if max_length < minimal_sv_length:
                        boundary_fails += 1
                        mutation["type"] = None
                        mutation["length"] = None
                        mutation["start"] = None
                        continue
                    
                    length_df = length_files.get(variant_type)
                    length = Selector.length(length_df, max_length, minimal_sv_length) if length_df is not None else minimal_sv_length
                    
                    if variant_type == "DUP":
                        for affected_node_id in node.get("affected_nodes", []):
                            affected_node = next((n for n in tree.get("nodes", []) if n["node"] == affected_node_id), None)
                            if affected_node:
                                affected_node["interval"][1] += length
                else:
                    mutation["type"] = None
                    mutation["length"] = None
                    mutation["start"] = None
                    continue
                
                mutation["length"] = length
                mutation["start"] = start_pos

        return tree, fails, boundary_fails
        
    except Exception as e:
        raise MSerror(f"Error augmenting mutations with length and position: {e}")


def process_single_tree(tree, variant_probabilities, minimal_sv_length, seed_offset, 
                        sv_length_files=None, tree_idx=None, assignments=None):
    """Processes a single tree."""
    random.seed(seed_offset)
    length_files = load_length_files_worker(sv_length_files)
    
    try:
        if assignments is not None and tree_idx is not None:
            augmented_tree = augment_type_from_assignments(tree, tree_idx, assignments)
        else:
            augmented_tree = augment_type(tree, variant_probabilities)
        
        augmented_tree, fails, boundary_fails = augment_length_and_position(
            augmented_tree, length_files, minimal_sv_length
        )
        return augmented_tree, fails, boundary_fails
    except MSerror as e:
        raise MSerror(f"Error processing the tree: {e}")


def process_batch(batch, variant_probabilities, num_threads, minimal_sv_length, 
                  seed_offset_start, sv_length_files=None, tree_idx_start=None, assignments=None):
    """Process a batch of trees in parallel."""
    if seed_offset_start is None:
        seed_list = [random.randint(0, 2**32 - 1) for _ in range(len(batch))]
    else:
        seed_list = [seed_offset_start + i for i in range(len(batch))]
    
    tree_indices = [tree_idx_start + i for i in range(len(batch))] if tree_idx_start is not None else [None] * len(batch)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(
            process_single_tree, 
            batch,
            [variant_probabilities] * len(batch), 
            [minimal_sv_length] * len(batch),
            seed_list,
            [sv_length_files] * len(batch),
            tree_indices,
            [assignments] * len(batch)
        ))
    return results


def main(json_file, output_json_file, demographic_file, chromosome, num_threads, batch_size=1000):
    """Processes trees in batches to handle millions of trees efficiently."""
    try:
        start_time = time.time()
        
        # Load all params from demographic file
        params = load_simulation_params(demographic_file)
        
        variant_probabilities = validate_sv_distribution(params['sv_distribution'])
        max_sv = validate_max_sv(params['max_sv'])
        sv_length_files = params['sv_length_files']
        minimal_sv_length = params['minimal_sv_length']
        readable_json = params['readable_json']
        seed = process_seed(params['seed'])
        
        MScompute(f"Generating variants for chromosome {chromosome} using {num_threads} threads")
        MScompute(f"Processing in batches of {batch_size} trees")
        
        sv_summary = ", ".join([f"{k}:{v}%" for k, v in variant_probabilities.items()])
        MScompute(f"SV distribution: {sv_summary}")
        
        if sv_length_files:
            MScompute(f"Custom length files: {sv_length_files}")

        if seed is not None:
            random.seed(seed)

        # Two-pass if max_sv is specified
        assignments = None
        if max_sv:
            MScompute("Using max_sv mode: two-pass processing")
            mutation_addresses, _ = collect_all_mutations(json_file)
            assignments, _ = assign_global_types(mutation_addresses, variant_probabilities, max_sv, seed)
            if seed is not None:
                random.seed(seed)

        # Initialize counters
        mutation_counts = {sv_type: 0 for sv_type in variant_probabilities.keys()}
        total_mutations = 0
        total_trees = 0
        total_fails = 0
        total_boundary_fails = 0
        
        with open(output_json_file, 'w') as outfile:
            outfile.write('[\n' if readable_json else '[')
            first_batch = True
            batch = []
            
            for tree in stream_json_trees(json_file):
                batch.append(tree)
                
                if len(batch) >= batch_size:
                    seed_offset = (seed + total_trees) if seed is not None else None
                    
                    results = process_batch(
                        batch, variant_probabilities, num_threads, 
                        minimal_sv_length, seed_offset, sv_length_files,
                        total_trees if assignments else None, assignments
                    )
                    
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
                    
                    MScompute(f"Processed {total_trees} trees...")
                    batch = []
            
            # Process remaining
            if batch:
                seed_offset = (seed + total_trees) if seed is not None else None
                results = process_batch(
                    batch, variant_probabilities, num_threads, 
                    minimal_sv_length, seed_offset, sv_length_files,
                    total_trees if assignments else None, assignments
                )
                
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
            
            outfile.write('\n]' if readable_json else ']')

        elapsed_time = time.time() - start_time
        
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
    parser.add_argument("--demographic", required=True, help="Path to the demographic JSON file with simulation params.")
    parser.add_argument("--output", required=True, help="Path to the output JSON file.")
    parser.add_argument("--chromosome", required=True, help="Chromosome name")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--batch_size", type=int, default=1000, help="Number of trees to process per batch")
    args = parser.parse_args()

    main(args.json, args.output, args.demographic, args.chromosome, args.threads, args.batch_size)