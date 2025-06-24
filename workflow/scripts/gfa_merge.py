"""
================================================================================
                    NODE MERGING OPTIMIZATION ALGORITHM
================================================================================

Author: Lucien Piat (INRAe, PangenOak Project)
Based on insights from: Sigfried Dubois

Description:
    This algorithm performs graph simplification by merging nodes that have 
    identical connectivity patterns. Two nodes can be merged when they form
    simple linear chains (e.g., A->B where A has only one outgoing edge to B
    and B has only one incoming edge from A).

Algorithm Complexity:
    - Original: O(n²) per iteration due to exhaustive pairwise checking
    - Optimized: O(n) per iteration using pre-computed connectivity patterns
    - Space complexity: O(n + e) where n = nodes, e = edges

Example:
    python gfa_merge.py input.gfa output.gfa --threads 4

Optimizations:
    - Pre-compute connectivity patterns to identify merge candidates in O(1)
    - Use multithreading for parallel processing of independent merge operations
    - Batch link updates to reduce overhead
    - Progress tracking with clear status updates

================================================================================
"""
#!/usr/bin/env python3
import argparse
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from threading import Lock
from typing import Dict, Set, Tuple, List
from io_handler import MSpangepopDataHandler, MScompute, MSsuccess, MSerror

@dataclass
class Segment:
    id: str
    sequence: str

@dataclass(frozen=True)
class Link:
    from_id: str
    from_orient: str
    to_id: str
    to_orient: str
    cigar: str

class GFAGraph:
    def __init__(self, num_threads=1):
        self.header = []
        self.segments = {}
        self.links = []
        self.walks = []
        self.other_lines = []
        self.num_threads = num_threads
        
        # For quick lookup - track incoming and outgoing connections separately
        self.outgoing = defaultdict(lambda: {"++": set(), "+-": set(), "-+": set(), "--": set()})
        self.incoming = defaultdict(lambda: {"++": set(), "+-": set(), "-+": set(), "--": set()})
        
        # Track which nodes have been merged into others
        self.merged_into = {}  # node_id -> node it was merged into
        
        # Thread safety
        self.merge_lock = Lock()
        
        # Optimization: Pre-compute merge candidates
        self.single_out_nodes = {"++": set(), "+-": set(), "-+": set(), "--": set()}
        self.single_in_nodes = {"++": set(), "+-": set(), "-+": set(), "--": set()}
        
    def parse_file(self, filename):
        """Parse GFA file"""
        print(f"Parsing GFA file")
        
        line_count = 0
        with open(filename, 'r') as f:
            for line in f:
                line_count += 1
                    
                if line.startswith('H'):
                    self.header.append(line.strip())
                elif line.startswith('S'):
                    parts = line.strip().split('\t')
                    seg_id = parts[1]
                    sequence = parts[2]
                    self.segments[seg_id] = Segment(seg_id, sequence)
                elif line.startswith('L'):
                    parts = line.strip().split('\t')
                    link = Link(parts[1], parts[2], parts[3], parts[4], parts[5])
                    self.links.append(link)
                elif line.startswith('W'):
                    self.walks.append(line.strip())
                else:
                    self.other_lines.append(line.strip())
        
        print(f"  Found {len(self.segments)} segments, {len(self.links)} links")
        
        # Build connection maps
        print("Building connection maps...")
        self._build_connections()
        
    def _build_connections(self):
        """Build connection maps and identify potential merge candidates"""
        
        for i, link in enumerate(self.links):
            if i % 100000 == 0 and i > 0:
                print(f"  Processing link {i}/{len(self.links)}...")
                
            orient_key = link.from_orient + link.to_orient
            self.outgoing[link.from_id][orient_key].add(link.to_id)
            
            # Build incoming connections
            if link.from_orient == '+' and link.to_orient == '+':
                self.incoming[link.to_id]['++'].add(link.from_id)
            elif link.from_orient == '+' and link.to_orient == '-':
                self.incoming[link.to_id]['+-'].add(link.from_id)
            elif link.from_orient == '-' and link.to_orient == '+':
                self.incoming[link.to_id]['-+'].add(link.from_id)
            else:  # -- 
                self.incoming[link.to_id]['--'].add(link.from_id)
        
        # Pre-compute nodes with single connections for faster merge detection
        self._update_single_connection_sets()
    
    def _update_single_connection_sets(self):
        """Update sets of nodes with single incoming/outgoing connections"""
        self.single_out_nodes = {"++": set(), "+-": set(), "-+": set(), "--": set()}
        self.single_in_nodes = {"++": set(), "+-": set(), "-+": set(), "--": set()}
        
        for node_id in self.segments:
            # Check outgoing connections
            out_conn = self.outgoing[node_id]
            total_out = sum(len(conn) for conn in out_conn.values())
            
            if total_out == 1:
                for orient, connections in out_conn.items():
                    if len(connections) == 1:
                        self.single_out_nodes[orient].add(node_id)
                        break
            
            # Check incoming connections
            in_conn = self.incoming[node_id]
            total_in = sum(len(conn) for conn in in_conn.values())
            
            if total_in == 1:
                for orient, connections in in_conn.items():
                    if len(connections) == 1:
                        self.single_in_nodes[orient].add(node_id)
                        break
    
    def find_mergeable_pairs(self):
        """Find pairs of nodes that can be merged based on the criteria - OPTIMIZED"""
        print("Finding mergeable pairs...")
        start_time = time.time()
        
        mergeable = []
        processed = set()
        
        # Only check nodes that have exactly one outgoing connection
        candidates = self.single_out_nodes['++'] | self.single_out_nodes['--']
        print(f"  Checking {len(candidates)} candidate nodes")
        
        checked = 0
        for node_id in candidates:
            if node_id in processed or node_id not in self.segments:
                continue
                
            checked += 1
            if checked % 100000 == 0:
                print(f"    Processed {checked}/{len(candidates)} candidates...")
            
            out_conn = self.outgoing[node_id]
            
            # Case 1: Node has exactly one ++ outgoing connection
            if len(out_conn['++']) == 1 and sum(len(conn) for orient, conn in out_conn.items() if orient != '++') == 0:
                neighbor_id = list(out_conn['++'])[0]
                if neighbor_id in processed or neighbor_id not in self.segments:
                    continue
                    
                # Quick check: is neighbor in single_in_nodes['++']?
                if neighbor_id in self.single_in_nodes['++']:
                    in_conn = self.incoming[neighbor_id]
                    if (len(in_conn['++']) == 1 and 
                        list(in_conn['++'])[0] == node_id and
                        sum(len(conn) for orient, conn in in_conn.items() if orient != '++') == 0):
                        
                        mergeable.append((node_id, neighbor_id, '++'))
                        processed.add(node_id)
                        processed.add(neighbor_id)
                        
            # Case 2: Node has exactly one -- outgoing connection  
            elif len(out_conn['--']) == 1 and sum(len(conn) for orient, conn in out_conn.items() if orient != '--') == 0:
                neighbor_id = list(out_conn['--'])[0]
                if neighbor_id in processed or neighbor_id not in self.segments:
                    continue
                    
                if neighbor_id in self.single_in_nodes['--']:
                    in_conn = self.incoming[neighbor_id]
                    if (len(in_conn['--']) == 1 and 
                        list(in_conn['--'])[0] == node_id and
                        sum(len(conn) for orient, conn in in_conn.items() if orient != '--') == 0):
                        
                        mergeable.append((node_id, neighbor_id, '--'))
                        processed.add(node_id)
                        processed.add(neighbor_id)
        
        find_time = time.time() - start_time
        print(f"Found {len(mergeable)} mergeable pairs")
        return mergeable
    
    def _merge_single_pair(self, pair_data):
        """Merge a single pair of nodes - thread-safe"""
        node1_id, node2_id, merge_type = pair_data
        
        with self.merge_lock:
            # Double-check nodes still exist (might have been merged by another thread)
            if node1_id not in self.segments or node2_id not in self.segments:
                return None
                
            seg1 = self.segments[node1_id]
            seg2 = self.segments[node2_id]
            
            # Create merged sequence
            if merge_type == '++':
                new_sequence = seg1.sequence + seg2.sequence
            elif merge_type == '--':
                new_sequence = seg1.sequence + self._reverse_complement(seg2.sequence)
            else:
                new_sequence = seg1.sequence + seg2.sequence
            
            # Update node1's sequence
            self.segments[node1_id].sequence = new_sequence
            
            # Remove node2
            del self.segments[node2_id]
            
            # Track the merge
            self.merged_into[node2_id] = node1_id
            
            return (node1_id, node2_id)
    
    def merge_nodes_parallel(self, mergeable_pairs):
        """Merge nodes in parallel"""
        if not mergeable_pairs:
            return 0
            
        print(f"Merging {len(mergeable_pairs)} pairs using {self.num_threads} threads...")
        start_time = time.time()
        
        merged_count = 0
        
        if self.num_threads == 1:
            # Single-threaded execution
            for i, pair in enumerate(mergeable_pairs):
                if i % 10000 == 0:
                    print(f"  Merged {i}/{len(mergeable_pairs)} pairs...")
                result = self._merge_single_pair(pair)
                if result:
                    merged_count += 1
        else:
            # Multi-threaded execution
            batch_size = max(1, len(mergeable_pairs) // (self.num_threads * 4))
            completed = 0
            
            with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
                # Submit batches
                futures = []
                for i in range(0, len(mergeable_pairs), batch_size):
                    batch = mergeable_pairs[i:i + batch_size]
                    for pair in batch:
                        futures.append(executor.submit(self._merge_single_pair, pair))
                
                # Process results
                for future in as_completed(futures):
                    result = future.result()
                    if result:
                        merged_count += 1
                    
                    completed += 1
                    if completed % 10000 == 0:
                        print(f"  Processed {completed}/{len(futures)} merge operations...")
        
        merge_time = time.time() - start_time
        print(f"  Merged {merged_count} pairs in {merge_time:.2f}s")
        
        return merged_count
    
    def update_links_batch(self):
        """Update all links after merging - optimized batch processing"""
        print("Updating links after merging...")
        start_time = time.time()
        
        new_links = []
        seen_links = set()
        
        for i, link in enumerate(self.links):
            if i % 100000 == 0 and i > 0:
                print(f"  Processing link {i}/{len(self.links)}...")
                
            new_link = None
            
            # Resolve merged nodes
            from_id = link.from_id
            to_id = link.to_id
            
            while from_id in self.merged_into:
                from_id = self.merged_into[from_id]
            while to_id in self.merged_into:
                to_id = self.merged_into[to_id]
            
            # Skip self-loops between merged nodes
            if from_id == to_id and link.from_id != link.to_id:
                continue
                
            new_link = Link(from_id, link.from_orient, to_id, link.to_orient, link.cigar)
            
            # Add the link if we haven't seen it before
            link_key = (new_link.from_id, new_link.from_orient, new_link.to_id, new_link.to_orient, new_link.cigar)
            if link_key not in seen_links:
                new_links.append(new_link)
                seen_links.add(link_key)
        
        self.links = new_links
        
        update_time = time.time() - start_time
        print(f"  Updated {len(new_links)} links in {update_time:.2f}s")
        
    def update_all_walks(self):
        """Update all walks after all merges are complete"""
        if not self.walks:
            return
            
        print(f"Updating {len(self.walks)} walks...")
        start_time = time.time()
        
        new_walks = []
        
        for i, walk in enumerate(self.walks):

            print(f"  Processing walk {i}/{len(self.walks)}...")
                
            parts = walk.split('\t')
            if len(parts) > 6:
                path = parts[6]
                
                # Parse path into nodes
                nodes = self._parse_path(path)
                
                # Update each node to its final merged version
                updated_nodes = []
                for node in nodes:
                    orientation = node[0]
                    node_id = node[1:]
                    
                    # Find the final node this was merged into
                    final_node_id = node_id
                    while final_node_id in self.merged_into:
                        final_node_id = self.merged_into[final_node_id]
                    
                    updated_nodes.append(orientation + final_node_id)
                
                # Collapse consecutive identical nodes
                collapsed_nodes = []
                prev = None
                for node in updated_nodes:
                    if node != prev:
                        collapsed_nodes.append(node)
                        prev = node
                
                # Reconstruct path
                parts[6] = ''.join(collapsed_nodes)
                
            new_walks.append('\t'.join(parts))
        
        self.walks = new_walks
        walk_time = time.time() - start_time
        print(f"  Updated walks in {walk_time:.2f}s")
    
    def _parse_path(self, path):
        """Parse a path string into a list of oriented nodes"""
        nodes = []
        i = 0
        while i < len(path):
            if path[i] in '><':
                orientation = path[i]
                i += 1
                node_id = ''
                while i < len(path) and path[i] not in '><':
                    node_id += path[i]
                    i += 1
                nodes.append(orientation + node_id)
            else:
                i += 1
        return nodes
    
    def _reverse_complement(self, seq):
        """Get reverse complement of a sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                     'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def rebuild_connections(self):
        """Rebuild connection maps after merging"""
        print("Rebuilding connection maps...")
        
        self.outgoing.clear()
        self.incoming.clear()
        
        for link in self.links:
            orient_key = link.from_orient + link.to_orient
            self.outgoing[link.from_id][orient_key].add(link.to_id)
            
            if link.from_orient == '+' and link.to_orient == '+':
                self.incoming[link.to_id]['++'].add(link.from_id)
            elif link.from_orient == '+' and link.to_orient == '-':
                self.incoming[link.to_id]['+-'].add(link.from_id)
            elif link.from_orient == '-' and link.to_orient == '+':
                self.incoming[link.to_id]['-+'].add(link.from_id)
            else:
                self.incoming[link.to_id]['--'].add(link.from_id)
        
        # Update single connection sets for next iteration
        self._update_single_connection_sets()
    
    def write_file(self, filename):
        """Write GFA file"""
        print(f"Writing output")
        
        with open(filename, 'w') as f:
            # Write headers
            for header in self.header:
                f.write(header + '\n')
            
            # Write segments (sorted numerically if possible)
            segments_written = 0
            for seg_id in sorted(self.segments.keys(), key=lambda x: int(x) if x.isdigit() else float('inf')):
                seg = self.segments[seg_id]
                f.write(f'S\t{seg.id}\t{seg.sequence}\n')
                segments_written += 1
                
            # Write links
            for i, link in enumerate(self.links):
                f.write(f'L\t{link.from_id}\t{link.from_orient}\t{link.to_id}\t{link.to_orient}\t{link.cigar}\n')
            
            # Write walks
            for walk in self.walks:
                f.write(walk + '\n')
            
            # Write other lines
            for line in self.other_lines:
                f.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge nodes in GFA file with multithreading support')
    parser.add_argument('input', help='Input GFA file')
    parser.add_argument('output', help='Output GFA file')
    parser.add_argument('--threads', type=int, default=1, 
                       help='Number of threads to use (default: 1)')
    args = parser.parse_args()
    
    if args.threads < 1:
        MSerror("Error: Number of threads must be at least 1")
        sys.exit(1)
    
    print("=" * 80)
    MScompute("Merging nodes, this can take a while")
    print("=" * 80)
    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Threads: {args.threads}")
    print()
    
    total_start_time = time.time()
    
    # Load and parse GFA file
    graph = GFAGraph(num_threads=args.threads)
    graph.parse_file(args.input)
    
    initial_nodes = len(graph.segments)
    initial_links = len(graph.links)
    print(f"\nInitial statistics:")
    print(f"  Nodes: {initial_nodes:,}")
    print(f"  Links: {initial_links:,}")
    
    total_merge_count = 0
    iteration = 1
    
    # Iterative merging until no more merges possible
    print(f"\n{'='*50}")
    MScompute("Pahse 1: Merging linear chains")
    print(f"{'='*50}")
    
    while True:
        print(f"\n--- Iteration {iteration} ---")
        iteration_start = time.time()
        
        # Find mergeable pairs
        mergeable_pairs = graph.find_mergeable_pairs()
        
        if not mergeable_pairs:
            print()
            MScompute("No more linear chains to merge.")
            break
        
        # Merge nodes
        merged_this_iteration = graph.merge_nodes_parallel(mergeable_pairs)
        total_merge_count += merged_this_iteration
        
        # Update links and rebuild connections
        graph.update_links_batch()
        graph.rebuild_connections()
        
        iteration_time = time.time() - iteration_start
        current_nodes = len(graph.segments)
        current_links = len(graph.links)
        
        print(f"  Iteration {iteration} summary:")
        print(f"    Nodes: {current_nodes:,} (reduced by {len(mergeable_pairs):,})")
        print(f"    Links: {current_links:,}")
        print(f"    Time: {iteration_time:.2f}s")
        
        iteration += 1
    
    # Update walks after all merges are complete
    print(f"\n{'='*50}")
    MScompute("Phase 2: Updating walks")
    print(f"{'='*50}")
    graph.update_all_walks()
    
    # Final statistics
    final_nodes = len(graph.segments)
    final_links = len(graph.links)
    reduction_nodes = (1 - final_nodes / initial_nodes) * 100 if initial_nodes > 0 else 0
    reduction_links = (1 - final_links / initial_links) * 100 if initial_links > 0 else 0
    
    total_time = time.time() - total_start_time
    
    print(f"\n{'='*50}")
    MScompute("Merge completed")
    print(f"{'='*50}")
    print(f"Initial nodes: {initial_nodes:,}")
    print(f"Final nodes: {final_nodes:,}")
    print(f"Nodes merged: {total_merge_count:,}")
    print(f"Node reduction: {reduction_nodes:.2f}%")
    print()
    print(f"Initial links: {initial_links:,}")
    print(f"Final links: {final_links:,}")
    print(f"Link reduction: {reduction_links:.2f}%")
    print()
    print(f"Total processing time: {total_time:.2f}s")
    print(f"Threads used: {args.threads}")
    
    graph.write_file(args.output)
    
    MSsuccess("Done!")
    print(f"\n{'='*80}")

if __name__ == "__main__":
    main()