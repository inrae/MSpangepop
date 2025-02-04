# Test script

# GFA format (simplified)
gfa_filename = os.path.join(output_dir, f"{chromosome_name}_msprime_simulation.gfa")
with open(gfa_filename, 'w') as f:
    #Header (why not)
    f.write("H\tVN:Z:1.0\n")

    #Print all nodes with no sequence
    for node in ts_chrom.nodes():
        f.write(f"S\t{node.id}\t*\n")
        
    #Add all edges
    for edge in ts_chrom.edges():
        f.write(f"L\t{edge.parent}\t+\t{edge.child}\t+\t0M\n")