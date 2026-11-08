# Going Even Deeper

## Readable JSON

mspangenome includes an option to output JSON files in human-readable format with proper indentation and formatting. While this makes the files larger and slower to process, it enables manual inspection and debugging.

## Markov Matrices for SNPs and Insertions

The graph code uses a substitution matrix to generate SNP bases and an order-1 Markov
matrix to generate insertion sequences. In `workflow/scripts/graph_classes.py`, the
currently selected SNP matrix is `random_matrix`, and the insertion matrix is
`simple_at_bias_matrix`, both defined in `workflow/scripts/matrix.py`. The alternative
`blosum_matrix` in that file is not currently selected by the graph code. You can change
the imports in `graph_classes.py` to select another matrix or edit the matrices to change
their substitution probabilities and sequence composition.

