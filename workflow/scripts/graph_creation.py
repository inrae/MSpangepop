"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

This script constructs a variation graph from genomic sequences and mutation data,
representing structural variations across multiple lineages in a unified graph structure.
The resulting graph captures SNPs, insertions, deletions, inversions, and duplications
while maintaining the relationships between different evolutionary lineages.

Workflow:
    The main() function orchestrates the following pipeline:
    
    1. INITIALIZATION PHASE:
       - Create tracking objects for mutations (MutationRecap) and variant 
         visualization (VariantSizeVisualizer)
       - Read input FASTA sequences (one per locus) and augmented traversal JSON
       - Calculate reference genome length from input sequences
       - Auto-scale parallelism based on available threads (workers = threads // 2)
    
    2. PARALLEL MEMORY-EFFICIENT GRAPH CONSTRUCTION (Phase 1):
       - Extract lineage information from traversal data
       - Process subgraphs in parallel using ProcessPoolExecutor
       - For each locus (in parallel): 
           * Initialize graph from sequence
           * Apply mutations (SNP, INS, DEL, INV, DUP)
           * Assign LOCAL node IDs (starting at 1 per subgraph)
           * Save to temp file as GFA fragment
           * Free graph from memory
       - Collect mutation/variant tracking data from workers
       - Calculate global ID offsets from node counts
    
    3. GRAPH MERGING WITH ID REMAPPING (Phase 2):
       - Stream through all temp files sequentially
       - Remap local node IDs to global IDs using cumulative offsets
       - Add connecting edges between consecutive subgraphs
       - Concatenate paths for each lineage across all subgraphs
       - Write final merged GFA file
    
    4. OUTPUT GENERATION (Phase 3):
       - Read final GFA and export lineage sequences as FASTA files
       - Apply reverse complement for nodes with "-" orientation
       - Generate comprehensive visualization plots
       - Write detailed mutation recap file with statistics
    
PARALLELIZATION:
    - Phase 1 is parallelized across subgraphs (significant speedup)
    - Phase 2 and 3 are sequential (I/O bound)
    - Thread scaling: --threads N → N/2 workers, each with N/(N/2) I/O threads

MEMORY MANAGEMENT:
    - Each subgraph is built, saved to temp file, then freed
    - Only metadata (SubgraphTempFile) kept in memory during Phase 1
    - Phase 2 streams from temp files, never loads all data at once
    
INVERSION HANDLING:
    - Inverted nodes are marked with "-" orientation in GFA paths
    - FASTA generation applies reverse complement for "-" oriented nodes

REQUIRED INPUTS:
    --splited_fasta       : Path to split FASTA file containing sequence segments
                           (one sequence per tree in the ARG)
    --augmented_traversal : Path to JSON file containing the ORDERED mutation list and lineages
                           information from ARG traversal and mutation augmentation
    --threads             : Total CPU threads available (auto-scales parallelism)
"""

import argparse
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List
from enum import Enum

from io_handler import MSpangepopDataHandler, MSerror, MSsuccess, MScompute
from graph_utils import (
    gather_lineages, MutationRecap, VariantSizeVisualizer, NodeIDAllocator
)
from graph_classes import Graph

class MutationType(Enum):
    """Enumeration of supported mutation types in the variation graph."""
    SNP = "SNP"
    INS = "INS"
    DEL = "DEL"
    INV = "INV"
    DUP = "DUP"

class SubgraphTempFile:
    """Tracks a saved subgraph temp file."""
    __slots__ = (
        "index", "path", "node_count",
        "local_start_id",   # Local ID of first node (for connecting edge)
        "local_end_id",     # Local ID of last node (for connecting edge)
        "mutations_data",   # Collected mutation records for merging into recap
        "variants_data"     # Collected variant records for merging into visualizer
    )
    
    def __init__(self, index: int, path: str, node_count: int,
                 local_start_id: int, local_end_id: int,
                 mutations_data: list, variants_data: list):
        self.index = index
        self.path = path
        self.node_count = node_count
        self.local_start_id = local_start_id
        self.local_end_id = local_end_id
        self.mutations_data = mutations_data
        self.variants_data = variants_data

# ============================================================================
# SUBGRAPH PROCESSING FUNCTIONS
# ============================================================================

def process_and_save_single_subgraph(args: tuple) -> SubgraphTempFile:
    """
    Process one subgraph: init -> apply mutations -> save to temp.
    Uses LOCAL node IDs (remapped to global during merge).
    
    Designed to run in a separate process (all args must be picklable).
    """
    (
        seq_data,          # (seq_id, seq_str) tuple - picklable
        tree_data,
        tree_index,
        lineages,          # set of lineage IDs
        temp_path,
        subgraph_index,
        sample,
        chromosome,
        io_threads
    ) = args
    
    # Reconstruct SeqRecord (Bio.SeqRecord isn't picklable)
    from Bio.Seq import Seq # type: ignore
    from Bio.SeqRecord import SeqRecord # type: ignore
    seq_record = SeqRecord(Seq(seq_data[1]), id=seq_data[0])
    
    # Local tracking (will be merged into main recap/visualizer later)
    local_mutations = []
    local_variants = []
    
    # 1. Initialize graph
    sequence = str(seq_record.seq)
    graph = Graph()
    graph.build_from_sequence(sequence, lineages)
    
    # 2. Apply mutations
    tree_interval = tree_data.get("initial_tree_interval", [0, 0])
    tree_start = tree_interval[0]
    
    for node_data in tree_data.get("nodes", []):
        node_id = node_data.get("node")
        mutations = node_data.get("mutations", [])
        affected_nodes = set(node_data.get("affected_nodes", []))
        
        if not mutations:
            continue
        
        tree_lineages_set = set(tree_data.get("lineages", []))
        affected_lineages = affected_nodes.intersection(tree_lineages_set)
        
        if not affected_lineages:
            continue
        
        for mutation in mutations:
            _apply_mutation_to_graph(
                graph=graph,
                mutation=mutation,
                tree_start=tree_start,
                tree_index=tree_index,
                node_id=node_id,
                affected_lineages=affected_lineages,
                local_mutations=local_mutations,
                local_variants=local_variants
            )
    
    # 4. Assign LOCAL node IDs (starting at 1 for each subgraph)
    local_allocator = NodeIDAllocator(start_id=1)
    assigned_nodes = set()
    
    local_start_id = None
    local_end_id = None
    
    for lineage, path in graph.paths.items():
        if lineage == "Ancestral":
            continue
        for edge in path.path_edges:
            if edge.node1 not in assigned_nodes:
                allocated_id = local_allocator.allocate(edge.node1)
                assigned_nodes.add(edge.node1)
                if local_start_id is None:
                    local_start_id = allocated_id
            if edge.node2 not in assigned_nodes:
                allocated_id = local_allocator.allocate(edge.node2)
                assigned_nodes.add(edge.node2)
                local_end_id = allocated_id
    
    # Assign IDs to any remaining nodes
    for node in graph.nodes:
        if node not in assigned_nodes:
            local_allocator.allocate(node)
            assigned_nodes.add(node)
    
    # Ensure we have valid start/end IDs
    if graph.start_node and graph.start_node.id:
        local_start_id = graph.start_node.id
    if graph.end_node and graph.end_node.id:
        local_end_id = graph.end_node.id
    
    # Default if still None
    if local_start_id is None:
        local_start_id = 1
    if local_end_id is None:
        local_end_id = max(1, local_allocator.current_id - 1)
    
    # 5. Save to temp file with local IDs
    MSpangepopDataHandler.save_subgraph_temp(graph, temp_path, sample, chromosome, max_workers=io_threads)
    
    # 6. Get node count before freeing
    node_count = len(graph.nodes)
    del graph
    
    return SubgraphTempFile(
        index=subgraph_index,
        path=temp_path,
        node_count=node_count,
        local_start_id=local_start_id,
        local_end_id=local_end_id,
        mutations_data=local_mutations,
        variants_data=local_variants
    )


def _apply_mutation_to_graph(
    graph: Graph,
    mutation: dict,
    tree_start: int,
    tree_index: int,
    node_id: int,
    affected_lineages: set,
    local_mutations: list,
    local_variants: list
):
    """Apply a single mutation with validation and tracking (parallel-safe)."""
    mut_type_str = mutation.get("type")
    try:
        mut_type = MutationType(mut_type_str) if mut_type_str is not None else None
    except ValueError:
        mut_type = None
    start = mutation.get("start")
    length = mutation.get("length")
    
    # Helper to record mutation
    def record_mutation(success: bool, error_msg: str = None, lineages: set = None):
        local_mutations.append({
            "tree_index": tree_index,
            "node_id": node_id,
            "type": mut_type.value if mut_type is not None else None,
            "position": start,
            "length": length,
            "affected_lineages": sorted(list(lineages or affected_lineages)),
            "success": success,
            "error": error_msg
        })
    
    # Helper to record variant
    def record_variant(success: bool, lineages: set = None):
        local_variants.append({
            "type": mut_type.value if mut_type is not None else None,
            "position": start,
            "length": length if mut_type != MutationType.SNP else 1,
            "success": success,
            "affected_lineages": list(lineages or affected_lineages)
        })
    
    # Handle None/skipped
    if mut_type is None:
        record_mutation(False, "Mutation skipped in previous step")
        return
    
    if start is None:
        record_mutation(False, "Mutation has no start position")
        record_variant(False)
        return
    
    relative_start = start - tree_start
    
    # Validate per lineage
    failed_lineages = []
    error_messages = []
    
    for lineage in affected_lineages:
        path = graph.paths.get(lineage)
        if not path:
            failed_lineages.append(lineage)
            error_messages.append(f"Lineage {lineage} not found")
            continue
        
        current_path_length = path.node_count + 1
        err = _validate_mutation(mut_type, relative_start, length, current_path_length)
        if err:
            failed_lineages.append(lineage)
            error_messages.append(err)
    
    valid_lineages = affected_lineages - set(failed_lineages)
    
    if failed_lineages:
        record_mutation(False, "; ".join(error_messages), set(failed_lineages))
        if not valid_lineages:
            record_variant(False)
            return
    
    # Apply mutation
    success = False
    error_msg = None
    
    try:
        if mut_type == MutationType.SNP:
            graph.add_snp(relative_start, valid_lineages)
        elif mut_type == MutationType.INS:
            graph.add_ins(relative_start, length, valid_lineages)
        elif mut_type == MutationType.DEL:
            graph.add_del(relative_start, relative_start + length, valid_lineages)
        elif mut_type == MutationType.INV:
            graph.add_inv(relative_start, relative_start + length - 1, valid_lineages)
        elif mut_type == MutationType.DUP:
            graph.add_tdup(relative_start, relative_start + length - 1, valid_lineages)
        else:
            raise MSerror(f"Unknown mutation type: {mut_type}")
        success = True
    except Exception as e:
        error_msg = str(e)
    
    record_mutation(success, error_msg, valid_lineages)
    record_variant(success, valid_lineages)


def _validate_mutation(mut_type: MutationType, relative_start: int, length: int, path_length: int):
    """Returns error string or None if valid."""
    if mut_type == MutationType.SNP:
        if relative_start == 0:
            return "SNP at first position would break connectivity"
        if relative_start >= path_length - 1:
            return f"SNP at/beyond last position (pos {relative_start}, len {path_length})"
    elif mut_type == MutationType.INS:
        if relative_start == 0:
            return "INS at position 0 would break connectivity"
        if relative_start >= path_length:
            return f"INS position {relative_start} >= path length {path_length}"
    elif mut_type == MutationType.DEL:
        if length is None:
            return "DEL has no length"
        if relative_start == 0:
            return "DEL at first position would break connectivity"
        if relative_start + length > path_length - 1:
            return "DEL would affect last position"
        if length >= path_length - 2:
            return "DEL would leave less than 3 nodes"
    elif mut_type in (MutationType.INV, MutationType.DUP):
        if length is None:
            return f"{mut_type.value} has no length"
        if relative_start == 0:
            return f"{mut_type.value} at first position would break connectivity"
        if relative_start + length - 1 >= path_length - 1:
            return f"{mut_type.value} would affect last position"
    return None

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main(splited_fasta: str, augmented_traversal: str, output_file: str,
         sample: str, chromosome: str, fasta_folder: str, tmp_folder: str, 
         recap_file: str = None, variant_plot_dir: str = None, 
         threads: int = 1) -> None:
    """
    Main function with automatic parallel processing based on available threads.
    """
    # Auto-scale: use half threads as workers, rest for I/O per worker
    # Minimum 1 worker, minimum 1 I/O thread per worker
    num_workers = max(1, threads // 2)
    io_threads_per_worker = max(1, threads // num_workers)
    
    # Initialize tracking
    recap = MutationRecap(sample, chromosome)
    
    # Read inputs
    MScompute("Reading input files")
    sequences = MSpangepopDataHandler.read_fasta(splited_fasta)
    traversal = MSpangepopDataHandler.read_json(augmented_traversal)
    
    reference_length = sum(len(str(record.seq)) for record in sequences)
    var_visualizer = VariantSizeVisualizer(sample, chromosome, reference_length)
    
    # Create temp directory
    temp_dir = os.path.join(tmp_folder, f"pangraph_{sample}_{chromosome}")
    os.makedirs(temp_dir, exist_ok=True)
    MScompute(f"Temp directory: {temp_dir}")
    
    try:
        tree_lineages = gather_lineages(traversal)
        
        if len(sequences) != len(tree_lineages):
            raise MSerror(f"Mismatch: {len(sequences)} sequences vs {len(tree_lineages)} trees")
        
        num_subgraphs = len(sequences)
        
        # Adjust workers if more workers than subgraphs
        actual_workers = min(num_workers, num_subgraphs)
        
        MScompute(f"Configuration: {actual_workers} workers, {io_threads_per_worker} I/O threads each")
        
        # Prepare picklable arguments
        tasks = []
        for i, (seq_record, (tree_index, lineages)) in enumerate(zip(sequences, tree_lineages)):
            seq_data = (seq_record.id, str(seq_record.seq))
            temp_path = os.path.join(temp_dir, f"subgraph_{i:06d}.gfa.tmp")
            
            tasks.append((
                seq_data,
                traversal[i],
                tree_index,
                lineages,
                temp_path,
                i,
                sample,
                chromosome,
                io_threads_per_worker
            ))
        
        # =====================================================================
        # PHASE 1: Process subgraphs in parallel with progress tracking
        # =====================================================================
        MScompute(f"Phase 1: Processing {num_subgraphs} subgraphs")
        
        temp_files: List[SubgraphTempFile] = [None] * len(tasks)
        completed = 0
        total = len(tasks)
        last_pct = -1
        
        MScompute(f"  Progress: 0/{total} (0%)")
        
        with ProcessPoolExecutor(max_workers=actual_workers) as executor:
            futures = {executor.submit(process_and_save_single_subgraph, task): task[5] for task in tasks}
            
            for future in as_completed(futures):
                try:
                    result = future.result()
                    temp_files[result.index] = result
                    completed += 1
                    
                    # Progress update every 2% or at completion
                    pct = (completed * 100) // total
                    if pct >= last_pct + 2 or completed == total:
                        last_pct = pct
                        MScompute(f"  Progress: {completed}/{total} ({pct}%)")
                        
                except Exception as e:
                    idx = futures[future]
                    raise MSerror(f"Subgraph {idx} failed: {e}")
        
        # =====================================================================
        # Merge tracking data from workers into main objects
        # =====================================================================
        MScompute("Merging tracking data from workers")
        
        for tf in temp_files:
            for mut_record in tf.mutations_data:
                recap.add_mutation(
                    tree_index=mut_record["tree_index"],
                    node_id=mut_record["node_id"],
                    mutation_type=mut_record["type"],
                    position=mut_record["position"],
                    length=mut_record["length"],
                    affected_lineages=set(mut_record["affected_lineages"]),
                    success=mut_record["success"],
                    error_msg=mut_record["error"]
                )
            
            for var_record in tf.variants_data:
                var_visualizer.add_variant(
                    variant_type=var_record["type"],
                    position=var_record["position"],
                    length=var_record["length"],
                    success=var_record["success"],
                    affected_lineages=set(var_record["affected_lineages"])
                )
            
            # Free memory
            tf.mutations_data = []
            tf.variants_data = []
        
        # Calculate global ID offsets
        id_offsets = [0]
        cumulative = 0
        for tf in temp_files:
            cumulative += tf.node_count
            id_offsets.append(cumulative)
        
        total_nodes = cumulative
        MSsuccess(f"Phase 1 complete: {len(temp_files)} subgraphs, {total_nodes:,} total nodes")
        
        # =====================================================================
        # PHASE 2: Merge temp files into final GFA with ID remapping
        # =====================================================================
        MScompute("Phase 2: Merging subgraphs with ID remapping")
        
        lineage_lengths = MSpangepopDataHandler.merge_temp_files_to_gfa(
            temp_files=temp_files,
            id_offsets=id_offsets,
            output_path=output_file,
            sample=sample,
            chromosome=chromosome
        )
        
        # =====================================================================
        # PHASE 3: Write FASTA from final GFA
        # =====================================================================
        MScompute("Phase 3: Writing FASTA")
        
        MSpangepopDataHandler.write_fasta_from_gfa(
            gfa_path=output_file,
            sample=sample,
            chromosome=chromosome,
            fasta_folder=fasta_folder
        )
        
        # Set lineage lengths for plots
        var_visualizer.set_lineage_lengths(lineage_lengths)
        
        MSsuccess("Graph construction complete!")
        
    except Exception as e:
        recap.summary["fatal_error"] = str(e)
        raise
    
    finally:
        # Cleanup temp files
        MScompute("Cleaning up temp files")
        shutil.rmtree(temp_dir, ignore_errors=True)
        
        # Save recap
        if recap_file:
            recap.save_recap(recap_file)
            MSsuccess("Saved recap")
        
        # Generate plots
        if variant_plot_dir:
            os.makedirs(variant_plot_dir, exist_ok=True)
            
            MScompute("Generating visualization plots")
            
            var_visualizer.save_size_distribution_plot(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_variant_sizes.png"))
            var_visualizer.save_variant_density_plot(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_variant_density.png"))
            var_visualizer.save_shared_variants_heatmap(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_shared_variants.png"))
            var_visualizer.save_variant_type_proportions_plot(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_variant_proportions.png"))
            var_visualizer.save_cumulative_variants_plot(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_cumulative_variants.png"))
            var_visualizer.save_lineage_lengths_plot(
                os.path.join(variant_plot_dir, f"{sample}_chr{chromosome}_graph_lineage_lengths.png"))
            
            MSsuccess("Visualization plots saved")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create variation graph from mutations")
    parser.add_argument("--splited_fasta", required=True, help="Path to split FASTA file")
    parser.add_argument("--augmented_traversal", required=True, help="Path to augmented traversal JSON")
    parser.add_argument("--output_file", required=True, help="Path to output GFA file")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--chromosome", required=True, help="Chromosome identifier")
    parser.add_argument("--fasta_folder", required=True, help="Directory to save all fasta")
    parser.add_argument("--tmp_folder", required=True, help="Directory for temporary files")
    parser.add_argument("--recap_file", help="Path to save recap file")
    parser.add_argument("--variant_plot_dir", help="Directory to save variant size plots")
    parser.add_argument("--threads", type=int, default=4, 
                        help="Total CPU threads available - auto-scales parallelism (default: 4)")

    args = parser.parse_args()

    main(
        splited_fasta=args.splited_fasta,
        augmented_traversal=args.augmented_traversal,
        output_file=args.output_file,
        sample=args.sample,
        chromosome=args.chromosome,
        fasta_folder=args.fasta_folder,
        tmp_folder=args.tmp_folder,
        recap_file=args.recap_file,
        variant_plot_dir=args.variant_plot_dir,
        threads=args.threads
    )