# Configuring the Pipeline

## Simulation Configuration

### Master Configuration File

The simulation parameters are defined in the `.config/masterconfig.yaml` file. This file controls all aspects of the msprime simulation for each sample.

#### Configuration Structure

The configuration file uses YAML format with a `samples` section where each sample is defined by a unique name and associated parameters.

```yaml
samples:        # <- Leave samples indent as is 
  sample_2:     # <- The first indent is you simulation sample
    params: A   # <- The second indent is the parameters used in the simulation
  sample_1:
    params: B   # <- The parameters are NOT shared between simulations 
```
 
⚠️ **WARNING**: Please do a test run with small population_size, mutation_rate, and recombination_rate first.

#### Sample Parameters

Each sample in the configuration requires the following parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `fasta_gz` | string | Path to the FASTA file containing the genome sequence |
| `chr_n` | integer | Number of chromosomes (contigs) in the FASTA file |
| `population_size` | integer | Effective population size (Ne) for the simulation |
| `mutation_rate` | float | Mutation rate (µ) - probability of mutations per base per generation |
| `recombination_rate` | float | Recombination rate (r) - frequency of recombination events |
| `sample_size` | integer | Number of individual samples to simulate |
| `model` | string | Mutation model for msprime: `'binary'`, `'infinite_alleles'`, or `'jc69'` |
| `minimal_variant_size` | integer | Minimum size for generated variants |
| `readable_json` | boolean | If `True`, outputs human-readable JSON (larger file size) |
| `seed` | integer/string/None | Random seed for reproducible simulations |

#### Global Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `memory_multiplier` | float | Multiplier for scaling memory allocation across all rules |
| `output_dir` | string | Directory where results will be stored |
| `sv_type_file` | string | Path to the structural variant type configuration file |

### Structural Variant Configuration

The distribution of structural variant types is specified in `.config/sv_type.yaml`. This file defines the percentage of each variant type in the simulation:

- **SNP**: Single Nucleotide Polymorphisms
- **DEL**: Deletions
- **INS**: Insertions
- **INV**: Inversions
- **DUP**: Duplications

⚠️ **Important**: The sum of all percentages must equal exactly 100.

### Variant Size Distribution

The size distribution for each structural variant type is defined in TSV files located in the `simulation_data/` directory:

- `size_distribDEL.tsv` - Deletion size distribution
- `size_distribDUP.tsv` - Duplication size distribution
- `size_distribINS.tsv` - Insertion size distribution
- `size_distribINV.tsv` - Inversion size distribution

Each TSV file should contain size bins and their corresponding probabilities for the respective variant type.

## Cluster Configuration

### SLURM Profile Configuration

The cluster execution parameters are defined in `.config/snakemake/profiles/slurm/config.yaml`. This file controls how jobs are submitted and executed on the SLURM cluster.

⚠️ **Note**: This file is intended for advanced users familiar with SLURM job scheduling. Modifying these settings can affect job execution, resource allocation, and queue behavior.

## Best Practices

1. **Test First**: Always run a small test simulation before launching large-scale analyses
2. **Parameter Scaling**: Start with conservative values for population size, mutation rate, and recombination rate
3. **Memory Management**: Adjust `memory_multiplier` if jobs are failing due to memory constraints

## Example config file 
```yaml
samples:
  s_aureus:
    fasta_gz: "staphylococcus_aureus_genome.fa.gz"
    chr_n: 1
    population_size: 5000
    mutation_rate: 0.5e-7     # <- Always start with a low value
    recombination_rate: 1e-10 # <- Always start with a low value
    sample_size: 5
    model: binary
    minimal_variant_size: 1
    readable_json: False
    seed: 1234
```