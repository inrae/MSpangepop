# Output Files

## Graph Files

The primary output of MSpangepop is a GFA (Graphical Fragment Assembly) file representing the variation graph. This file contains nodes representing DNA sequences and edges representing connections between them, along with walk paths showing how each simulated individual traverses the graph.

## Result folder outline : 
```bash
results/../
├── 01_plots                                 # <- Many plots
├── 02_simulation_data 
│   ├── chr_1        
│   │   ├── graph_creation_recap.txt         # <- Recap for the graph creation step              
│   │   ├── msprime_recap.txt                # <- Recap for the msprime simulation step         
│   │   ├── ancestry_ts.trees                # <- the .ts files contains the raw msprime simulation
│   │   └── mutated_ts.trees                           
│   └── small_test_global_recap.txt          # <- Global recap         
├── 03_graph    
│   ├── chr_1_graph_polished.gfa             # <- Final GFA                
│   └── traversal                           
│       └── chr_1_augmented_traversal.json   # <- The file used for graph creation    
└── benchmark                                                       
```