# Mutation log

Every `graph_creation` run writes one registry per chromosome:

    results/<sample>/simulation_data/chr_<N>/<sample>_chr_<N>_graph_creation_recap.tsv

It records every mutation the simulator attempted, in the order it was applied.

## Format

A TSV with metadata in leading `#`.

```python
import pandas as pd
mut = pd.read_csv(path, sep="\t", comment="#")
```

```r
mut <- read.delim(path, comment.char = "#")
```

| column | meaning |
|---|---|
| `order` | 1-based rank of the attempt, the order they were applied in |
| `tree_index` | index of the local genealogy, one per recombination block |
| `tree_node` | node of that genealogy carrying the mutation |
| `type` | `SNP`, `DEL`, `INS`, `INV` or `TDUP`, empty when the draw step skipped it |
| `position` | start, in chromosome coordinates |
| `locus_offset` | start, relative to the locus the graph engine works in |
| `length` | affected length in bases, 1 for a SNP |
| `lineages` | every affected lineage, comma separated |
| `status` | `OK` or `FAILED` |
| `error` | why it failed, empty otherwise |

Missing values are empty fields.

## Header lines

    #format       mspangenome-mutation-log/1
    #sample       SA_A_rep1
    #chromosome   1
    #generated    2026-10-31T10:19:55
    #duration     0:04:12
    #attempted    3026
    #successful   2945
    #failed       81
    #type         DEL   412   405   7

One `#type` line per mutation type, with attempted, successful and failed counts.
`#fatal_error` appears only when the run aborted.

## Order of application

Mutations are applied in a preorder walk of the genealogy, so a mutation on a node
is always applied before any mutation below it.

## Positions are not reference coordinates

`position` is where the mutation was drawn, in the coordinates of the sequence the
locus started from. Each applied mutation moves the ones that follow it on the same
lineage, so an insertion at offset 100 shifts every later position in the lineages
that carry it. This is what makes nested variants possible. It means the columns cannot be read as coordinates on a finished haplotype. To follow a locus through its mutations one at a time, see
[locus_replay](https://github.com/Lucien-Piat/MSpangenome_methods).
