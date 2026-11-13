# Going Even Deeper

## Readable JSON

MSpangenome includes an option to output JSON files in human-readable format with proper indentation and formatting. While this makes the files larger and slower to process, it enables manual inspection and debugging.

## Markov Matrices for SNPs and Insertions

The graph code uses a substitution matrix to generate SNP bases and an order-1 Markov
matrix to generate insertion sequences. In `workflow/scripts/graph_classes.py`, the
currently selected SNP matrix is `random_matrix`, and the insertion matrix is
`simple_at_bias_matrix`, both defined in `workflow/scripts/matrix.py`. The alternative
`blosum_matrix` in that file is not currently selected by the graph code. You can change
the imports in `graph_classes.py` to select another matrix or edit the matrices to change
their substitution probabilities and sequence composition.

## Graph Path Operations

MSpangenome represents variants by changing how selected lineage paths traverse the graph. An operation can target one or several paths, allowing variants to be unique to a lineage or shared across lineages. Edges record the sides used to exit and enter each node, preserving orientation (++, +-, -+, --). Operations can be combined or overlap to represent compound variants, such as a deletion within an inversion.

| Operation | Function | Used for | Path change |
|-----------|----------|----------|-------------|
| `bypass(a,b)` | Skip nodes `a` to `b` | Deletions | Creates a shortcut edge |
| `loop(a,b)` | Duplicate nodes `a` to `b` | Tandem duplications | Adds a loop-back and repeat |
| `invert(a,b)` | Reverse nodes `a` to `b` | Inversions | Reverses the path direction |
| `swap(a,node)` | Replace the node at position `a` | SNPs | Substitutes a single node |
| `paste(a,a+1,nodes)` | Insert between adjacent nodes | Insertions | Adds a new node sequence |

New nodes, such as those used for insertions, are generated with an order-1 Markov model.

