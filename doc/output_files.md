# Output Files

The main deliverable is the chromosome GFA in `03_graph/`. The workflow also writes
tree-sequence files, run recaps, optional plots, per-rule logs, and benchmark files.

## Output layout

With the default `output_dir: "results/"` and sample `sample_name`, the persistent output
tree is approximately:

```text
results/
└── sample_name/
    ├── 01_plots/
    │   └── chr_1/
    │       ├── 01_tree_plots/              # unless succinct: true
    │       ├── 02_arg_plots/               # unless succinct: true
    │       ├── 03_graph_plots/             # unless succinct: true
    │       └── visualization_complete.txt  # omitted in succinct mode
    ├── 02_simulation_data/
    │   ├── sample_name_global_recap.txt
    │   └── chr_1/
    │       ├── sample_name_chr_1_ancestry_ts.trees
    │       ├── sample_name_chr_1_full_ancestry_ts.trees
    │       ├── sample_name_chr_1_mutated_ts.trees
    │       ├── sample_name_chr_1_msprime_recap.txt
    │       └── sample_name_chr_1_graph_creation_recap.txt
    ├── 03_graph/
    │   └── chr_1/
    │       ├── sample_name_chr_1_graph.gfa
    │       ├── sample_name_chr_1_graph_stats.txt
    │       ├── fasta/
    │       │   └── sample_name_chr1.fasta.gz
    │       └── traversal/
    │           └── sample_name_chr_1_augmented_traversal.json
    └── benchmark/
        └── chr_1/                          # one benchmark file per rule
logs/                                       # default log_dir, outside results/
└── sample_name/
    ├── setup.log
    └── chr_1/                              # one log file per rule
```

The actual base directory for rule logs is the configured `log_dir` (default `logs/`),
which is separate from `output_dir`. The corresponding paths are
`logs/sample_name/setup.log` and `logs/sample_name/chr_1/<rule>.log`.

If `replicates` is greater than one, `sample_name` is suffixed with `_rep1`, `_rep2`, etc.
The sample output folder and its logs use that expanded name.

## Files and purpose

| File or directory | Purpose | Persistence |
|---|---|---|
| `03_graph/chr_<n>/<sample>_chr_<n>_graph.gfa` | Final variation graph for a chromosome | Persistent |
| `03_graph/chr_<n>/<sample>_chr_<n>_graph_stats.txt` | Node and edge counts before and after final unchopping | Persistent |
| `03_graph/chr_<n>/fasta/` | FASTA sequences reconstructed from graph paths | Persistent |
| `03_graph/chr_<n>/traversal/` | Variant-annotated traversal JSON used to create the graph | Persistent |
| `02_simulation_data/chr_<n>/*.trees` | Mutated, ancestry, and full-ARG tree sequences from msprime/tskit | Persistent |
| `02_simulation_data/chr_<n>/*recap.txt` | Per-chromosome simulation and graph-construction summaries | Persistent |
| `02_simulation_data/<sample>_global_recap.txt` | Configuration, software and run information | Persistent |
| `01_plots/` | Tree, ARG, and graph visualizations | Omitted when `succinct: true` |
| `<log_dir>/<sample>/` | Standard output and error for each workflow rule | Persistent |
| `benchmark/chr_<n>/` | Snakemake timing/resource benchmark files | Persistent |

The exact plot files depend on the visualizer. The GFA subgraphs created in the temporary
directory are Snakemake temporary outputs and are removed after the merge rule completes.
Other intermediate files are also marked temporary in the workflow. The temporary directory
name defaults to `00_temp/` and can be changed with `tempfile_location`.

The range-expansion step writes its expanded YAML to the configured `expanded_config` path
(default `.config/expanded_config.yaml`) and its per-sample demographic JSON files under
`.config/expanded_demographics/`. These configuration files are outside `output_dir`.
