# Input Files

## Primary Input

The input file of the workflow is a telomere to telomere assembly of the studied organism. 

The wokflow will run a simulation on each contig in the fasta file so a fragmented assembly will produce a lot of GFA file. 

The alphabet of the assembly must be {A, T, C, G}, no Ns. 

This file should be compressed with bgzip for optimal performance.
