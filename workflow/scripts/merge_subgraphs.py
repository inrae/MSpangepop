#!/usr/bin/env python3
"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Concatenate the locus subgraphs into the final chromosome graph.

Loci are independent by construction: `split_recombination` cuts the
reference at recombination breakpoints and locus boundaries are kept free of
mutations, so a subgraph can be compacted without looking at its neighbours.
Each `graph_creation` worker therefore unchops its own locus through
`run_unchop` below, as soon as it has built it, and throws the chopped copy
away. vg only ever works on a few kilobases, and the chopped form exists on
disk only for the loci currently in flight.

Steps:
    1. Concatenate the unchopped subgraphs: shift node IDs into a common space,
       join consecutive loci with an edge, and glue each lineage path back
       together.
    2. Unchop once more to merge the locus junctions.
"""

import argparse
import glob
import os
import shlex
import subprocess

from io_handler import MSerror, MSsuccess, MScompute

# vg holds the graph in PackedGraph form (-p) and writes paths as P-lines (-W),
# which is what the concatenation below parses.
UNCHOP_PIPELINE = (
    "set -o pipefail; "
    "vg convert -g {input} -t {threads} -p | "
    "vg mod -u --remove-non-path -t {threads} - | "
    "vg convert -fW -t {threads} - > {output}"
)


def run_unchop(input_path: str, output_path: str, threads: int) -> None:
    """Run one vg unchop pipeline, raising with vg's own message on failure."""
    command = UNCHOP_PIPELINE.format(
        input=shlex.quote(input_path),
        output=shlex.quote(output_path),
        threads=threads
    )
    result = subprocess.run(["bash", "-c", command], capture_output=True, text=True)
    if result.returncode != 0:
        raise MSerror(
            f"Unchop of {os.path.basename(input_path)} failed: {result.stderr.strip()}"
        )


def concatenate_subgraphs(fragments: list, output_path: str) -> tuple:
    """
    Concatenate unchopped subgraphs into a single chromosome graph.

    Node IDs are local to each subgraph, so every ID is shifted by the number of
    IDs the previous subgraphs consumed. Consecutive loci are joined by an edge
    from the last node of one locus to the first node of the next, and each
    lineage path is the concatenation of its per-locus paths.

    Returns (node_count, edge_count).
    """
    paths_by_lineage = {}
    connecting_edges = []
    expected_paths = None
    offset = 0
    previous_end = None
    node_count = 0
    edge_count = 0

    with open(output_path, "w", buffering=8 * 1024 * 1024) as out_f:
        out_f.write("H\tVN:Z:1.1\n")

        for index, fragment in enumerate(fragments):
            max_local_id = 0
            first_segment = None
            last_segment = None
            paths_here = set()

            with open(fragment, "r", buffering=4 * 1024 * 1024) as in_f:
                for line in in_f:
                    if line.startswith("S\t"):
                        parts = line.split("\t", 2)
                        local_id = int(parts[1])
                        max_local_id = max(max_local_id, local_id)
                        out_f.write(f"S\t{local_id + offset}\t{parts[2]}")
                        node_count += 1

                    elif line.startswith("L\t"):
                        parts = line.rstrip("\n").split("\t")
                        node1, node2 = int(parts[1]), int(parts[3])
                        max_local_id = max(max_local_id, node1, node2)
                        out_f.write(
                            f"L\t{node1 + offset}\t{parts[2]}\t"
                            f"{node2 + offset}\t{parts[4]}\t{parts[5]}\n"
                        )
                        edge_count += 1

                    elif line.startswith("P\t"):
                        parts = line.rstrip("\n").split("\t")
                        name = parts[1]
                        segments = parts[2].split(",")
                        max_local_id = max(
                            max_local_id, max(int(s[:-1]) for s in segments)
                        )
                        shifted = [f"{int(s[:-1]) + offset}{s[-1]}" for s in segments]

                        # Locus boundaries carry no mutation, so every lineage
                        # enters and leaves a locus through the same two nodes.
                        # If that no longer holds the chromosome path would be
                        # stitched through the wrong node.
                        if first_segment is None:
                            first_segment, last_segment = shifted[0], shifted[-1]
                        elif (shifted[0], shifted[-1]) != (first_segment, last_segment):
                            raise MSerror(
                                f"Lineages disagree on the boundaries of {fragment}: "
                                f"{name} runs {shifted[0]}..{shifted[-1]} instead of "
                                f"{first_segment}..{last_segment}"
                            )

                        paths_here.add(name)
                        paths_by_lineage.setdefault(name, []).extend(shifted)

            if first_segment is None:
                raise MSerror(f"Subgraph {fragment} holds no path")

            if expected_paths is None:
                expected_paths = paths_here
            elif paths_here != expected_paths:
                raise MSerror(
                    f"Subgraph {fragment} does not carry the same lineages as the "
                    f"first one: missing {sorted(expected_paths - paths_here)}, "
                    f"unexpected {sorted(paths_here - expected_paths)}"
                )

            if previous_end is not None:
                connecting_edges.append(
                    f"L\t{previous_end[:-1]}\t{previous_end[-1]}\t"
                    f"{first_segment[:-1]}\t{first_segment[-1]}\t0M\n"
                )
            previous_end = last_segment
            offset += max_local_id

            if (index + 1) % 100 == 0 or index == len(fragments) - 1:
                MScompute(f"  Concatenated {index + 1}/{len(fragments)} subgraphs")

        out_f.write("".join(connecting_edges))
        edge_count += len(connecting_edges)

        for name, segments in paths_by_lineage.items():
            out_f.write(f"P\t{name}\t{','.join(segments)}\t*\n")

    return node_count, edge_count


def count_gfa(gfa_path: str) -> tuple:
    """Count nodes and edges of a GFA, the two figures `vg stats -z` reports."""
    nodes = 0
    edges = 0
    with open(gfa_path, "r", buffering=4 * 1024 * 1024) as f:
        for line in f:
            if line.startswith("S\t"):
                nodes += 1
            elif line.startswith("L\t"):
                edges += 1
    return nodes, edges


def main(subgraph_dir: str, output_file: str, stats_file: str,
         tmp_folder: str, threads: int) -> None:
    fragments = sorted(glob.glob(os.path.join(subgraph_dir, "subgraph_*.gfa")))
    if not fragments:
        raise MSerror(f"No subgraph found in {subgraph_dir}")

    concatenated = os.path.join(tmp_folder, "concatenated_graph.gfa")
    os.makedirs(tmp_folder, exist_ok=True)

    try:
        MScompute(f"Concatenating {len(fragments)} subgraphs")
        merged_nodes, merged_edges = concatenate_subgraphs(fragments, concatenated)
        MSsuccess(f"Concatenated graph: {merged_nodes:,} nodes, {merged_edges:,} edges")

        # The per-locus passes could not touch the junctions between loci, this
        # one does. It runs on the compacted graph, not on the chopped one.
        MScompute("Unchopping locus junctions")
        run_unchop(concatenated, output_file, threads)
        final_nodes, final_edges = count_gfa(output_file)

        with open(stats_file, "w") as f:
            chopped_record = os.path.join(subgraph_dir, "chopped_nodes.txt")
            if os.path.exists(chopped_record):
                with open(chopped_record) as chopped:
                    f.write("[gfa_merge] Initial GFA state :\n")
                    f.write(f"nodes\t{chopped.read().strip()}\n")
            f.write(f"[gfa_merge] Subgraphs :\n{len(fragments)}\n")
            f.write("[gfa_merge] Concatenated, before junction unchop :\n")
            f.write(f"nodes\t{merged_nodes}\nedges\t{merged_edges}\n")
            f.write("[gfa_merge] Final state :\n")
            f.write(f"nodes\t{final_nodes}\nedges\t{final_edges}\n")

        MSsuccess(f"Final GFA: {final_nodes:,} nodes, {final_edges:,} edges")

    finally:
        if os.path.exists(concatenated):
            os.remove(concatenated)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate the locus subgraphs and unchop their junctions")
    parser.add_argument("--subgraph_dir", required=True,
                        help="Directory holding the chopped subgraph GFA files")
    parser.add_argument("--output_file", required=True, help="Path to the final GFA")
    parser.add_argument("--stats_file", required=True, help="Path to the stats file")
    parser.add_argument("--tmp_folder", required=True,
                        help="Directory for intermediate files")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for the junction unchop (default: 1)")

    args = parser.parse_args()

    main(
        subgraph_dir=args.subgraph_dir,
        output_file=args.output_file,
        stats_file=args.stats_file,
        tmp_folder=args.tmp_folder,
        threads=args.threads
    )
