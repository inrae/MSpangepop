"""
================================================================================
                    NODE MERGING OPTIMIZATION ALGORITHM
================================================================================

Author: Lucien Piat (INRAe, PangenOak Project)
Based on insights from: Sigfried Dubois

Description:
    This algorithm performs graph simplification by merging nodes that have 
    identical connectivity patterns. Two nodes can be merged when:

Algorithm Complexity:

Example:

================================================================================
"""
#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Set, Tuple, List

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
    def __init__(self):
        self.header = []
        self.segments = {}
        self.links = []
        self.walks = []
        self.other_lines = []
        
        # For quick lookup - track incoming and outgoing connections separately
        self.outgoing = defaultdict(lambda: {"++": set(), "+-": set(), "-+": set(), "--": set()})
        self.incoming = defaultdict(lambda: {"++": set(), "+-": set(), "-+": set(), "--": set()})
        
        # Track which nodes have been merged into others
        self.merged_into = {}  # node_id -> node it was merged into
        
    def parse_file(self, filename):
        """Parse GFA file"""
        with open(filename, 'r') as f:
            for line in f:
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
        
        # Build connection maps
        for link in self.links:
            orient_key = link.from_orient + link.to_orient
            self.outgoing[link.from_id][orient_key].add(link.to_id)
            
            # Reverse orientation for incoming
            if link.from_orient == '+' and link.to_orient == '+':
                self.incoming[link.to_id]['++'].add(link.from_id)
            elif link.from_orient == '+' and link.to_orient == '-':
                self.incoming[link.to_id]['+-'].add(link.from_id)
            elif link.from_orient == '-' and link.to_orient == '+':
                self.incoming[link.to_id]['-+'].add(link.from_id)
            else:  # -- 
                self.incoming[link.to_id]['--'].add(link.from_id)
    
    def find_mergeable_pairs(self):
        """Find pairs of nodes that can be merged based on the criteria"""
        mergeable = []
        processed = set()
        
        total_nodes = len(self.segments)
        checked = 0
        
        for node_id in sorted(self.segments.keys()):
            if node_id in processed:
                continue
                
            checked += 1
            if checked % 100 == 0:
                print(f"Checking node {checked}/{total_nodes} for merge candidates...")
            
            # Get all outgoing connections from this node
            out_conn = self.outgoing[node_id]
            
            # Case 1: Node has exactly one ++ outgoing connection
            if (len(out_conn['++']) == 1 and 
                len(out_conn['+-']) == 0 and 
                len(out_conn['-+']) == 0 and 
                len(out_conn['--']) == 0):
                
                neighbor_id = list(out_conn['++'])[0]
                if neighbor_id in processed:
                    continue
                    
                # Check if neighbor has exactly one ++ incoming connection (from this node)
                in_conn = self.incoming[neighbor_id]
                if (len(in_conn['++']) == 1 and 
                    list(in_conn['++'])[0] == node_id and
                    len(in_conn['+-']) == 0 and 
                    len(in_conn['-+']) == 0 and 
                    len(in_conn['--']) == 0):
                    
                    # This forms a simple linear path that can be merged
                    mergeable.append((node_id, neighbor_id, '++'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
                    
            # Case 2: Node has exactly one -- outgoing connection  
            elif (len(out_conn['--']) == 1 and 
                  len(out_conn['++']) == 0 and 
                  len(out_conn['+-']) == 0 and 
                  len(out_conn['-+']) == 0):
                
                neighbor_id = list(out_conn['--'])[0]
                if neighbor_id in processed:
                    continue
                    
                in_conn = self.incoming[neighbor_id]
                if (len(in_conn['--']) == 1 and 
                    list(in_conn['--'])[0] == node_id and
                    len(in_conn['++']) == 0 and 
                    len(in_conn['+-']) == 0 and 
                    len(in_conn['-+']) == 0):
                    
                    mergeable.append((node_id, neighbor_id, '--'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
        
        return mergeable
    
    def find_bidirectional_pairs(self):
        """Find pairs of nodes with mutual inverted links"""
        mergeable = []
        processed = set()
        
        for node_id in sorted(self.segments.keys()):
            if node_id in processed:
                continue
                
            out_conn = self.outgoing[node_id]
            
            # Check for bidirectional connections
            # Case 1: A->B (++) and B->A (--)
            for neighbor_id in out_conn['++']:
                if neighbor_id in processed:
                    continue
                neighbor_out = self.outgoing[neighbor_id]
                if node_id in neighbor_out['--']:
                    # Found bidirectional pair - merge them regardless of other connections
                    mergeable.append((node_id, neighbor_id, 'bidirectional'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
                    break
            
            # Case 2: A->B (--) and B->A (++)
            for neighbor_id in out_conn['--']:
                if neighbor_id in processed:
                    continue
                neighbor_out = self.outgoing[neighbor_id]
                if node_id in neighbor_out['++']:
                    # Found bidirectional pair - merge them regardless of other connections
                    mergeable.append((node_id, neighbor_id, 'bidirectional'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
                    break
            
            # Case 3: A->B (+-) and B->A (-+)
            for neighbor_id in out_conn['+-']:
                if neighbor_id in processed:
                    continue
                neighbor_out = self.outgoing[neighbor_id]
                if node_id in neighbor_out['-+']:
                    mergeable.append((node_id, neighbor_id, 'bidirectional'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
                    break
            
            # Case 4: A->B (-+) and B->A (+-)
            for neighbor_id in out_conn['-+']:
                if neighbor_id in processed:
                    continue
                neighbor_out = self.outgoing[neighbor_id]
                if node_id in neighbor_out['+-']:
                    mergeable.append((node_id, neighbor_id, 'bidirectional'))
                    processed.add(node_id)
                    processed.add(neighbor_id)
                    break
        
        return mergeable
    
    def merge_nodes(self, node1_id, node2_id, merge_type):
        """Merge two nodes into one"""
        if node1_id not in self.segments or node2_id not in self.segments:
            return None
            
        seg1 = self.segments[node1_id]
        seg2 = self.segments[node2_id]
        
        # Create merged segment
        if merge_type == '++':
            new_sequence = seg1.sequence + seg2.sequence
        elif merge_type == '--':
            new_sequence = seg1.sequence + self._reverse_complement(seg2.sequence)
        elif merge_type == 'bidirectional':
            # For bidirectional, just concatenate (or you could choose a different strategy)
            new_sequence = seg1.sequence + seg2.sequence
        else:  # '++/--'
            new_sequence = seg1.sequence + seg2.sequence
        
        # Update node1's sequence
        self.segments[node1_id].sequence = new_sequence
        
        # Remove node2
        del self.segments[node2_id]
        
        # Track the merge
        self.merged_into[node2_id] = node1_id
        
        # Update links
        new_links = []
        for link in self.links:
            # Skip links between node1 and node2
            if (link.from_id == node1_id and link.to_id == node2_id) or \
               (link.from_id == node2_id and link.to_id == node1_id):
                continue
                
            # Update links from node2 to come from node1
            elif link.from_id == node2_id:
                new_link = Link(node1_id, link.from_orient, link.to_id, link.to_orient, link.cigar)
                new_links.append(new_link)
                
            # Update links to node2 to go to node1
            elif link.to_id == node2_id:
                new_link = Link(link.from_id, link.from_orient, node1_id, link.to_orient, link.cigar)
                new_links.append(new_link)
            else:
                new_links.append(link)
        
        self.links = new_links
        
        return node1_id
    
    def update_all_walks(self):
        """Update all walks after all merges are complete"""
        new_walks = []
        
        for walk in self.walks:
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
    
    def write_file(self, filename):
        """Write GFA file"""
        with open(filename, 'w') as f:
            for header in self.header:
                f.write(header + '\n')
            
            for seg_id in sorted(self.segments.keys(), key=lambda x: int(x) if x.isdigit() else float('inf')):
                seg = self.segments[seg_id]
                f.write(f'S\t{seg.id}\t{seg.sequence}\n')
            
            for link in self.links:
                f.write(f'L\t{link.from_id}\t{link.from_orient}\t{link.to_id}\t{link.to_orient}\t{link.cigar}\n')
            
            for walk in self.walks:
                f.write(walk + '\n')
            
            for line in self.other_lines:
                f.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge nodes in GFA file')
    parser.add_argument('input', help='Input GFA file')
    parser.add_argument('output', help='Output GFA file')
    args = parser.parse_args()
    
    print("Loading GFA file...")
    graph = GFAGraph()
    graph.parse_file(args.input)
    
    initial_nodes = len(graph.segments)
    print(f"Initial number of nodes: {initial_nodes}")
    
    merge_count = 0
    iteration = 1
    
    # First pass: merge linear chains
    print("\n=== Phase 1: Merging linear chains ===")
    while True:
        print(f"\nIteration {iteration}:")
        print("Finding mergeable pairs...")
        mergeable_pairs = graph.find_mergeable_pairs()
        
        if not mergeable_pairs:
            print("No more linear chains to merge.")
            break
        
        print(f"Found {len(mergeable_pairs)} mergeable pairs")
        
        for i, (node1, node2, merge_type) in enumerate(mergeable_pairs):
            if i % 20 == 0:
                print(f"Merging pair {i+1}/{len(mergeable_pairs)}...")
            graph.merge_nodes(node1, node2, merge_type)
            merge_count += 1
        
        graph.rebuild_connections()
        iteration += 1
    
    # First pass: merge linear chains
    print("\n=== Phase 2: Merging bidirectional pairs ===")
    while True:
        print(f"\nIteration {iteration}:")
        print("Finding bidirectional pairs...")
    
        bidirectional_pairs = graph.find_bidirectional_pairs()
        
        if bidirectional_pairs:
            print(f"Found {len(bidirectional_pairs)} bidirectional pairs")
            for i, (node1, node2, merge_type) in enumerate(bidirectional_pairs):
                print(f"Merging bidirectional pair: {node1} <-> {node2}")
                graph.merge_nodes(node1, node2, merge_type)
                merge_count += 1
            graph.rebuild_connections()
            iteration += 1
        else:
            print("No more pairs found.")
            break
            

    # Update all walks after all merges are complete
    print("\nUpdating walks...")
    graph.update_all_walks()
    
    final_nodes = len(graph.segments)
    reduction = (1 - final_nodes / initial_nodes) * 100
    
    print(f"\nMerge complete!")
    print(f"Initial nodes: {initial_nodes}")
    print(f"Final nodes: {final_nodes}")
    print(f"Nodes merged: {merge_count}")
    print(f"Node reduction: {reduction:.2f}%")
    
    print(f"\nWriting output to {args.output}...")
    graph.write_file(args.output)
    print("Done!")

if __name__ == "__main__":
    main()