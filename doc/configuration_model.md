#  Creating a Demographic Model JSON for MSpangenome

Population structure, demographic events, and mutation/recombination rates are translated
into msprime inputs. `generation_time` is recorded as metadata; simulation times and rates
are interpreted per generation. The `simulation_params` section contains settings for the
MSpangenome workflow itself. For details about the demographic concepts supported by
msprime, see the [msprime demography documentation](https://tskit.dev/msprime/docs/stable/demography.html).

If you don't want to delve into the **lengthy** MSprime documentation, please stick with the default **Panmictic** Model, adapt mutation rate, recombination_rate and sample_size with your values.



## 1. File Structure Overview

A complete demographic model contains:
- **Simulation parameters** (genome file, chromosomes, SV distribution)
- **Evolutionary parameters** (mutation rate, recombination rate)
- **Population structure** (populations, samples, events)

## Evolutionary Parameters for Common Model Organisms

| Organism | Mutation Rate (per bp per gen) | Recombination Rate (per bp per gen) | Key References |
|----------|-------------------------------|-------------------------------------|----------------|
| **Human** | 1.0–1.5 × 10⁻⁸ | ~1 × 10⁻⁸ | Kong et al. 2012 *Nature*; Narasimhan et al. 2017 *Nat Commun* |
| **Drosophila melanogaster** | 2–5 × 10⁻⁹ | ~2.5 × 10⁻⁸ | Wang et al. 2023 *Genome Res*; Comeron et al. 2012 *PLOS Genet* |
| **Oak (*Quercus*)** | ~1 × 10⁻⁸ | ~1 × 10⁻⁸ | Sork et al. 2022 *Nat Commun*; Plomion et al. 2018 *Nat Plants* |
| **Yeast (*S. cerevisiae*)** | 1–4 × 10⁻¹⁰ * | ~3 × 10⁻⁶ * | Lang & Murray 2008 *Genetics* |
| **E. coli** | ~5 × 10⁻¹⁰ * | N/A (rare) | Drake 1991 *PNAS*; Lee et al. 2012 *PNAS* |

###  Basic Template

```json
{
  "name": "My_Model",
  "description": "Description of this model",
  
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {
      "SNP": 50,
      "DEL": 20,
      "INS": 15,
      "INV": 10,
      "DUP": 5
    },
    "seed": "optional_seed",
    "minimal_sv_length": 1
  },
  
  "evolutionary_params": {
    "mutation_rate": 1e-7,
    "recombination_rate": 1e-8,
    "generation_time": 25
  },
  
  "populations": [
    {"id": "pop1", "initial_size": 5000}
  ],
  
  "samples": [
    {"population": "pop1", "sample_size": 10}
  ]
}
```

## 2. Using Parameter Ranges

Ranges are supported only for the fields processed by `sample_ranges.py`: values under
`evolutionary_params`; population `initial_size`; the supported demographic-event fields
(`time`, `size`, `growth_rate`, `proportion`, `rate`); migration-matrix `rate` values;
`simulation_params.sv_distribution`, `max_sv`, and `minimal_sv_length`. Other numeric
fields are not expanded automatically. Each replicate samples its own values for supported
ranges before the workflow starts.

### Fixed vs Range Parameters

```json
{
  "evolutionary_params": {
    "mutation_rate": 1e-7,
    "recombination_rate": {"min": 1e-9, "max": 1e-7}
  }
}
```

### Available Distributions

| Distribution | Keywords | Description |
|-------------|----------|-------------|
| **Uniform** | `"uniform"` | Equal probability across range |
| **Log-uniform** | `"log_uniform"`, `"log"` | Uniform in log space (good for rates) |
| **Normal** | `"normal"` | Bell curve, truncated to range |
| **Auto** | `"auto"` or omitted | Intended to choose log-uniform for rates < 0.01 and uniform otherwise; see the note below |

**Implementation note:** supported fields are currently passed to the sampler with
prefixed names such as `evol_mutation_rate` and `event_0_rate`, so the automatic rate-name
check does not recognize them. In the current workflow, `auto` therefore resolves to
uniform for supported ranges. Specify `"distribution": "log_uniform"` explicitly when
log-uniform sampling is required; its minimum must be greater than zero.

```json
{
  "evolutionary_params": {
    "mutation_rate": {
      "min": 1e-8,
      "max": 1e-6,
      "distribution": "log_uniform"
    }
  }
}
```

---

## 3. Simulation Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `fasta_gz` | string | Reference FASTA path (plain or gzip-compressed) | `"genome.fa.gz"` |
| `chr_n` | integer | Number of chromosomes to simulate | `1` |
| `sv_distribution` | dict | Variant-type proportions, including SNPs (normalized to 100%) | See below |
| `max_sv` | dict | **Maximum count per SV type per chromosome** | See below |
| `sv_length_files` | dict | **Custom length distribution TSV files** | See below |
| `minimal_sv_length` | integer | Minimum length for generated structural variants | `1` |
| `seed` | string/int | Seed used by some stochastic steps (not a complete workflow seed; see below) | `"myseed"` or `12345` |
| `readable_json` | boolean | Human-readable JSON output | `false` |

The model's `seed` is read by variant generation. In the standard range-expansion path,
the seed is not copied to the expanded sample configuration used by the msprime rule, and
the range-expansion step is not given a seed. It therefore does **not** currently make the
whole workflow reproducible.

### SV Distribution

Controls the **proportion** of each variant type, including SNPs (normalized to 100%):

```json
{
  "simulation_params": {
    "sv_distribution": {
      "SNP": 50,
      "DEL": 20,
      "INS": 15,
      "INV": 10,
      "DUP": 5
    }
  }
}
```

With ranges:
```json
{
  "simulation_params": {
    "sv_distribution": {
      "SNP": {"min": 40, "max": 60},
      "DEL": {"min": 15, "max": 25},
      "INS": {"min": 10, "max": 20},
      "INV": {"min": 5, "max": 15},
      "DUP": {"min": 2, "max": 8}
    }
  }
}
```

### Max SV (Per-Chromosome Limits)

Limits the **number per variant type on each chromosome**. Mutations from the chromosome's
trees are randomly selected up to these limits:

```json
{
  "simulation_params": {
    "max_sv": {
      "SNP": 100,
      "DEL": 20,
      "INS": 20,
      "INV": 20,
      "DUP": 20
    }
  }
}
```

With ranges:
```json
{
  "simulation_params": {
    "max_sv": {
      "SNP": {"min": 50, "max": 200},
      "DEL": {"min": 10, "max": 50},
      "INS": {"min": 10, "max": 50},
      "INV": {"min": 5, "max": 30},
      "DUP": {"min": 5, "max": 30}
    }
  }
}
```

**Note:** If `max_sv` is specified, excess mutations are removed. Selection is randomized across the trees for that chromosome.

**Note:** A maximum is a cap, not a target. If the simulation produces fewer eligible mutations, the cap will not be reached.
### SV Length Files

Custom TSV files for length distributions (defaults used if omitted):

```json
{
  "simulation_params": {
    "sv_length_files": {
      "DEL": "my_data/del_lengths.tsv",
      "INS": "my_data/ins_lengths.tsv",
      "INV": "my_data/inv_lengths.tsv",
      "DUP": "my_data/dup_lengths.tsv"
    }
  }
}
```

Default files are in `simulation_data/size_distrib{TYPE}.tsv`. Each TSV must have a
`size_interval` column (intervals such as `[100, 200]`) and a `pb` column containing the
probability for that interval. Paths are resolved relative to the workflow's working
directory.

---

## 4. Evolutionary Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `mutation_rate` | float | Per-base mutation rate per-generation|
| `recombination_rate` | float | Per-base recombination rate per-generation |
| `generation_time` | integer | Years per generation |

```json
{
  "evolutionary_params": {
    "mutation_rate": 1e-8,
    "recombination_rate": 1e-8,
    "generation_time": 25
  }
}
```

**Note:** Choose `mutation_rate` carefully: high values can produce many mutations and make simulations expensive (for example, `1e-6` may be high for some genomes).

**Note:** Inversions, duplications and deletions must fit within the current locus. High recombination rates create shorter loci, which can limit or truncate the lengths available for these variants.

---

## 5. Populations and Samples

Here you can specify each population you want to simulate. Keep in mind that all individuals must coalesce to a common ancestor at some point.

```json
{
  "populations": [
    {"id": "pop1", "initial_size": 5000},
    {"id": "pop2", "initial_size": {"min": 1000, "max": 3000}}
  ],
  "samples": [
    {"population": "pop1", "sample_size": 10},
    {"population": "pop2", "sample_size": 5}
  ]
}
```

An optional top-level `migration_matrix` can set migration rates at time zero. Each entry
contains `time`, `source`, `dest`, and `rate`; the current loader applies entries whose
`time` is `0`. For changes at later times, use an `add_migration_rate_change` event.
---

## 6. Demographic Events

Here you can add events that will affect your population

### Population Size Change
```json
{
  "demographic_events": [
    {
      "type": "population_parameters_change",
      "time": 1000,
      "population": "pop1",
      "size": {"min": 5000, "max": 15000}
    }
  ]
}
```

### Migration
```json
{
  "demographic_events": [
    {
      "type": "add_migration_rate_change",
      "time": 500,
      "source": "pop1",
      "dest": "pop2",
      "rate": {"min": 1e-5, "max": 1e-3}
    }
  ]
}
```

### Mass Migration
```json
{
  "demographic_events": [
    {
      "type": "mass_migration",
      "time": 1500,
      "source": "pop1",
      "dest": "pop2",
      "proportion": 0.2
    }
  ]
}
```

---

## 7. Complete Example

The repository contains additional models in `./simulation_data`. Adapt the genome path and parameters to your data.
```json
{
  "name": "Example_Model",
  
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {
      "SNP": 50,
      "DEL": 20,
      "INS": 15,
      "INV": 10,
      "DUP": 5
    },
    "max_sv": {
      "SNP": 100,
      "DEL": 20,
      "INS": 20,
      "INV": 20,
      "DUP": 20
    },
    "sv_length_files": {
      "DEL": "simulation_data/size_distribDEL.tsv",
      "INS": "simulation_data/size_distribINS.tsv",
      "INV": "simulation_data/size_distribINV.tsv",
      "DUP": "simulation_data/size_distribDUP.tsv"
    },
    "seed": 42
  },
  
  "evolutionary_params": {
    "mutation_rate": 1e-8,
    "recombination_rate": 1e-8,
    "generation_time": 25
  },
  
  "populations": [
    {"id": "pop1", "initial_size": 5000}
  ],
  
  "samples": [
    {"population": "pop1", "sample_size": 3}
  ]
}
```

The following population structure illustrates a split from `POP_A` into `POP_B`. Samples refer to the declared population IDs.
```json
{
  "name": "Two_Population_Split",
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {"SNP": 50, "DEL": 20, "INS": 15, "INV": 10, "DUP": 5}
  },
  "evolutionary_params": {
    "mutation_rate": 1e-8,
    "recombination_rate": 1e-8,
    "generation_time": 25
  },
  "populations": [
    {"id": "POP_A", "initial_size": 5000},
    {"id": "POP_B", "initial_size": 1000}
  ],
  "samples": [
    {"population": "POP_A", "sample_size": 8},
    {"population": "POP_B", "sample_size": 3}
  ],
  "migration_matrix": [],
  "demographic_events": [
    {
      "type": "mass_migration",
      "time": 3000,
      "source": "POP_B",
      "dest": "POP_A",
      "proportion": 1.0
    }
  ]
}
```
<img src="demo.png">

---

## ✅ 8. Validation Checklist

- ✔ All populations have `id` and `initial_size`
- ✔ Every sample references an existing population
- ✔ `simulation_params` includes `fasta_gz`, `chr_n`, and `sv_distribution`
- ✔ `evolutionary_params` includes `mutation_rate` and `recombination_rate`
- ✔ SV distributions sum to approximately 100% (will be normalized)
- ✔ All ranges have `min` ≤ `max`
- ✔ If using `max_sv`, values are non-negative integers (`0` disables a variant type)
- ✔ If using `sv_length_files`, paths exist and are valid TSV files
