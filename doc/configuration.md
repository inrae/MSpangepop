# Configuring the Pipeline

## Simulation Configuration

### Master Configuration File

The simulation parameters are defined in the `.config/masterconfig.yaml` file. This file controls all aspects of the simulation pipeline, including individual samples and parameter sweeps.

#### Configuration Structure

The configuration file uses YAML format with two main sections:
- `samples`: Individual samples with fixed parameters
- `sweep_samples`: Automated parameter space exploration

```yaml
# Individual samples with fixed parameters
samples:
  sample_1:     # Unique sample identifier
    params: A   # Parameters for this specific sample

# Parameter sweeps for automated exploration
sweep_samples:
  sweep_1:      # Unique sweep identifier
    params: B   # Base parameters and sweep ranges
```

⚠️ **WARNING**: Please do a test run with small population sizes and low mutation/recombination rates first.

### Individual Samples

Each sample in the `samples` section requires the following parameters:

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `fasta_gz` | string | Path to compressed FASTA file containing the genome sequence | Required |
| `chr_n` | integer | Number of chromosomes (contigs) to simulate from the FASTA file | Required |
| `demographic_file` | string | Path to demographic JSON file containing population parameters | Required |
| `sv_distribution` | dict | SV type percentages: `{SNP: %, DEL: %, INS: %, INV: %, DUP: %}` | Required |
| `model` | string | Mutation model: `'binary'`, `'infinite_alleles'`, or `'jc69'` | `'binary'` |
| `minimal_sv_length` | integer | Minimum length for structural variants in base pairs | `1` |
| `readable_json` | boolean | If `True`, outputs human-readable JSON (larger file size) | `False` |
| `seed` | int/string/None | Random seed for reproducible simulations | `None` |

#### Example Individual Sample
```yaml
samples:
  s_aureus_baseline:
    fasta_gz: "staphylococcus_aureus_genome.fa.gz"
    chr_n: 1
    demographic_file: "simulation_data/Panmictic_Model.json"
    sv_distribution: {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0}
    seed: 1234
```

### Parameter Sweeps

The `sweep_samples` section enables automated exploration of parameter space by generating multiple samples with different parameter combinations.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `base_demographic_file` | string | Base demographic JSON file to modify | Required |
| `fasta_gz` | string | Path to compressed FASTA file | Required |
| `chr_n` | integer | Number of chromosomes to simulate | Required |
| `sweep_params` | dict | Demographic parameters to sweep (see below) | `{}` |
| `sv_sweep` | list | List of SV distributions to test | Base distribution |
| `sv_distribution` | dict | Default SV distribution if `sv_sweep` not specified | Required if no `sv_sweep` |
| `replicates` | integer | Number of replicates per parameter combination | `1` |
| `seed` | int/string/None | Fixed seed for ALL expanded samples | `None` |
| `model` | string | Mutation model | `'binary'` |
| `minimal_sv_length` | integer | Minimum SV length | `1` |
| `readable_json` | boolean | Human-readable JSON output | `False` |

#### Sweep Parameters

The `sweep_params` dictionary can include any parameter from the demographic JSON file:

- `mutation_rate`: List of mutation rates to test (e.g., `[1e-6, 1e-5, 1e-4]`)
- `recombination_rate`: List of recombination rates to test
- Any other demographic parameter from the JSON file

#### Example Parameter Sweep
```yaml
sweep_samples:
  mutation_sv_exploration:
    base_demographic_file: "simulation_data/Panmictic_Model.json"
    fasta_gz: "small_test_genome.fa.gz"
    chr_n: 1
    sweep_params:
      mutation_rate: [1e-7, 1e-6, 1e-5]
      recombination_rate: [1e-8, 1e-7]
    sv_sweep:
      - {SNP: 90, DEL: 5, INS: 5, INV: 0, DUP: 0}   # SNP-heavy
      - {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0} # Balanced
      - {SNP: 20, DEL: 30, INS: 30, INV: 15, DUP: 5} # SV-heavy
    replicates: 3
    seed: "fixed_seed"  # Same seed for all combinations
```

This example creates 3 × 2 × 3 × 3 = 54 simulations (3 mutation rates × 2 recombination rates × 3 SV distributions × 3 replicates).

### SV Type Distribution

Structural variant distributions are defined directly in the configuration as dictionaries:

```yaml
sv_distribution: {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0}
```

- **SNP**: Single Nucleotide Polymorphisms (%)
- **DEL**: Deletions (%)
- **INS**: Insertions (%)
- **INV**: Inversions (%)
- **DUP**: Duplications (%)

⚠️ **Important**: The sum of all percentages must equal 100. The pipeline will warn and normalize if not.

### Variant Size Distribution

The size distribution for each structural variant type is defined in TSV files located in the `simulation_data/` directory:

- `size_distribDEL.tsv` - Deletion size distribution
- `size_distribDUP.tsv` - Duplication size distribution
- `size_distribINS.tsv` - Insertion size distribution
- `size_distribINV.tsv` - Inversion size distribution

Each TSV file contains size bins and their corresponding probabilities.

### Global Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `memory_multiplier` | float | Memory scaling factor for all rules | `1` |
| `output_dir` | string | Base directory for simulation outputs | `"results/"` |



## Running with Parameter Sweeps

When parameter sweeps are defined, the pipeline automatically:

1. **Expands configurations**: Creates all parameter combinations
2. **Generates demographic files**: Modified versions saved to `.config/expanded_demographics/`
3. **Names samples systematically**: `sweepname_param1value_param2value_rep#`
4. **Handles seeds properly**: 
   - If seed specified: ALL samples use the SAME seed (isolates parameter effects)
   - If no seed: Each run gets a different random seed (true replicates)


## Cluster Configuration

### SLURM Profile Configuration

The cluster execution parameters are defined in `.config/snakemake/profiles/slurm/config.yaml`. This file controls how jobs are submitted and executed on the SLURM cluster.

⚠️ **Note**: This file is intended for advanced users familiar with SLURM job scheduling.

## Best Practices

1. **Test First**: Always run small test simulations before large-scale analyses
2. **Parameter Scaling**: Start with conservative mutation/recombination rates
3. **Memory Management**: Adjust `memory_multiplier` if jobs fail due to memory constraints
4. **Sweep Design**: 
   - Use fixed seeds to isolate parameter effects
   - Include replicates to assess variability
   - Start with coarse parameter grids, then refine
5. **SV Distributions**: Ensure percentages sum to 100 to avoid normalization warnings

## Complete Example Configuration

```yaml
# .config/masterconfig.yaml

# Individual samples
samples:
  baseline_simulation:
    fasta_gz: "small_test_genome.fa.gz"
    chr_n: 1
    demographic_file: "simulation_data/Island_Model.json"
    sv_distribution: {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0}
    seed: "baseline"

# Parameter sweeps
sweep_samples:
  comprehensive_sweep:
    base_demographic_file: "simulation_data/Panmictic_Model.json"
    fasta_gz: "small_test_genome.fa.gz"
    chr_n: 1
    sweep_params:
      mutation_rate: [1e-7, 1e-6, 1e-5]
    sv_sweep:
      - {SNP: 100, DEL: 0, INS: 0, INV: 0, DUP: 0}  # Pure SNPs
      - {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0} # Mixed
    replicates: 5
    seed: "sweep_seed"  # Fixed seed for all

# Global settings
memory_multiplier: 1
output_dir: "results/"
tempfile_location: "00_temp/"
```
