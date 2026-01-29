#  Creating a Demographic Model JSON for MSpangepop

**Warning:** Parameters defined in this config file are parsed and fed to MSprime. Describing all of them **exceeds** the scope of this documentation so please refer **directly** to MSprime documentation for further reference: [https://tskit.dev/msprime/docs/stable/demography.html](https://tskit.dev/msprime/docs/stable/demography.html)

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

Any numeric parameter can be specified as a range. Each replicate will randomly sample values.

### Fixed vs Range Parameters

```json
{
  "evolutionary_params": {
    "mutation_rate": 1e-7,                              // Fixed value
    "recombination_rate": {"min": 1e-9, "max": 1e-7}   // Range
  }
}
```

### Available Distributions

| Distribution | Keywords | Description |
|-------------|----------|-------------|
| **Uniform** | `"uniform"` | Equal probability across range |
| **Log-uniform** | `"log_uniform"`, `"log"` | Uniform in log space (good for rates) |
| **Normal** | `"normal"` | Bell curve, truncated to range |
| **Auto** | `"auto"` or omitted | Log-uniform for rates < 0.01, else uniform |

```json
"mutation_rate": {
  "min": 1e-8,
  "max": 1e-6,
  "distribution": "log_uniform"
}
```

---

## 3. Simulation Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `fasta_gz` | string | Compressed reference genome | `"genome.fa.gz"` |
| `chr_n` | integer | Number of chromosomes to simulate | `1` |
| `sv_distribution` | dict | Structural variant percentages (normalized to 100%) | See below |
| `max_sv` | dict | **Maximum count per SV type (global limit)** | See below |
| `sv_length_files` | dict | **Custom length distribution TSV files** | See below |
| `seed` | string/int | Random seed (optional) | `"myseed"` or `12345` |
| `readable_json` | boolean | Human-readable JSON output | `false` |

### SV Distribution

Controls the **proportion** of each variant type (normalized to 100%):

```json
"sv_distribution": {
  "SNP": 50,
  "DEL": 20,
  "INS": 15,
  "INV": 10,
  "DUP": 5
}
```

With ranges:
```json
"sv_distribution": {
  "SNP": {"min": 40, "max": 60},
  "DEL": {"min": 15, "max": 25},
  "INS": {"min": 10, "max": 20},
  "INV": {"min": 5, "max": 15},
  "DUP": {"min": 2, "max": 8}
}
```

### Max SV (Global Limits)

Limits the **total number** of each variant type across all trees. Mutations are randomly selected up to these limits:

```json
"max_sv": {
  "SNP": 100,
  "DEL": 20,
  "INS": 20,
  "INV": 20,
  "DUP": 20
}
```

With ranges:
```json
"max_sv": {
  "SNP": {"min": 50, "max": 200},
  "DEL": {"min": 10, "max": 50},
  "INS": {"min": 10, "max": 50},
  "INV": {"min": 5, "max": 30},
  "DUP": {"min": 5, "max": 30}
}
```

**Note:** If `max_sv` is specified, mutations exceeding the limits are deleted. Selection is randomized across all trees to avoid bias.

**Note:** If `max_sv` is specified, mutations and `mutation_rate` is low, the maximum number of SVs may not be reached.
### SV Length Files

Custom TSV files for length distributions (defaults used if omitted):

```json
"sv_length_files": {
  "DEL": "my_data/del_lengths.tsv",
  "INS": "my_data/ins_lengths.tsv",
  "INV": "my_data/inv_lengths.tsv",
  "DUP": "my_data/dup_lengths.tsv"
}
```

Default files are in `simulation_data/size_distrib{TYPE}.tsv`. 

---

## 4. Evolutionary Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `mutation_rate` | float | Per-base mutation rate per-generation|
| `recombination_rate` | float | Per-base recombination rate per-generation |
| `generation_time` | integer | Years per generation |

```json
"evolutionary_params": {
  "mutation_rate": 1e-8,
  "recombination_rate": 1e-8,
  "generation_time": 25
}
```

**Note:** Be carefull with the `mutation_rate` parameter, the simulation may be realy heavy if this number is to high (1e-6 for example)

**Note:** SV length for inversions, duplications and deletions cannot exceed the the current locus size, a high `recombination_rate` will truncate most SV. Lower the `recombination_rate` if you want large variations to appear.  

---

## 5. Populations and Samples

Here you can specify each population you want to simulate. Keep in mind that all individuals must coalesce to a common ancestor at some point.

```json
"populations": [
  {"id": "pop1", "initial_size": 5000},
  {"id": "pop2", "initial_size": {"min": 1000, "max": 3000}}
],

"samples": [
  {"population": "pop1", "sample_size": 10},
  {"population": "pop2", "sample_size": 5}
]
```
---

## 6. Demographic Events

Here you can add events that will affect your population

### Population Size Change
```json
{
  "type": "population_parameters_change",
  "time": 1000,
  "population": "pop1",
  "size": {"min": 5000, "max": 15000}
}
```

### Migration
```json
{
  "type": "add_migration_rate_change",
  "time": 500,
  "source": "pop1",
  "dest": "pop2",
  "rate": {"min": 1e-5, "max": 1e-3}
}
```

### Mass Migration
```json
{
  "type": "mass_migration",
  "time": 1500,
  "source": "pop1",
  "dest": "pop2",
  "proportion": 0.2
}
```

---

## 7. Complete Example

Below a complete model, you can find additional exampkles in the `./simulation_data` folder. You can try them on your genome and 
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
    "mutation_rate": 1e-6,
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

Graphical example :   
```json
"populations": [
    { "id": "POP_A", "initial_size": 5000 },
    { "id": "POP_B", "initial_size": 1000 }
  ],

  "samples": [
    { "population": "POP_1", "sample_size": 8 },
    { "population": "POP_2", "sample_size": 3 }
  ],
  
  "migration_matrix": [],

  "demographic_events": [
    {
      "type": "mass_migration",
      "time": 3000,
      "source": "POP_B",
      "dest": "POP_A",
      "proportion": 1.0,
      "_comment": "Population split (backward in time: POP_B merges into POP_A)"
    }
  ]
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
- ✔ If using `max_sv`, values are positive integers
- ✔ If using `sv_length_files`, paths exist and are valid TSV files