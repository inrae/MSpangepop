# Creating a Demographic Model JSON for MSpangepop


## 📝 1. File Structure Overview

A complete demographic model contains:
- **Simulation parameters** (genome file, chromosomes, SV distribution)
- **Evolutionary parameters** (mutation rate, recombination rate)
- **Population structure** (populations, samples, events)

### Basic Template

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
    "model": "binary",
    "minimal_sv_length": 1,
    "readable_json": false
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


## 🎲 2. Using Parameter Ranges

Any numeric parameter can be either fixed or specified as a range. When using ranges, each replicate will randomly sample values.

### Fixed vs Range Parameters

```json
{
  "evolutionary_params": {
    "mutation_rate": 1e-7,                        // Fixed value
    "recombination_rate": {"min": 1e-9, "max": 1e-7},  // Range
    "generation_time": 25                         // Fixed value
  }
}
```

### Range Syntax

Use `{"min": value, "max": value}` for any numeric parameter:

```json
{
  "populations": [
    {
      "id": "pop1",
      "initial_size": {"min": 1000, "max": 10000}  // Population size range
    }
  ],
  
  "demographic_events": [
    {
      "type": "add_migration_rate_change",
      "time": 1000,
      "source": "pop1",
      "dest": "pop2",
      "rate": {"min": 1e-5, "max": 1e-3}  // Migration rate range
    }
  ]
}
```

### SV Distribution Ranges

Each SV type can have a range. The distribution is normalized to 100% after sampling:

```json
"sv_distribution": {
  "SNP": {"min": 30, "max": 50},
  "DEL": {"min": 15, "max": 25},
  "INS": {"min": 10, "max": 20},
  "INV": {"min": 10, "max": 20},
  "DUP": {"min": 5, "max": 15}
}
```

---

## 3. Parameter Sections

### Simulation Parameters

Located in `"simulation_params"`:

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `fasta_gz` | string | Compressed reference genome | `"genome.fa.gz"` |
| `chr_n` | integer | Number of chromosomes to simulate (cant be greater than the number fo contigs present in `fasta_gz`)| `1` |
| `sv_distribution` | dict | Structural variant percentages | See above |
| `seed` | string/int | Random seed (optional) | `"myseed"` or `12345` |
| `readable_json` | boolean | Human-readable output | `false` |

### Evolutionary Parameters

Located in `"evolutionary_params"`:

| Parameter | Type | Description | Range Example |
|-----------|------|-------------|---------------|
| `mutation_rate` | float | Per-base mutation rate | `{"min": 1e-8, "max": 1e-6}` |
| `recombination_rate` | float | Per-base recombination rate | `{"min": 1e-9, "max": 1e-7}` |
| `generation_time` | integer | Years per generation | `25` |

---

## 4. Demographic Events

### Population Size Change

```json
{
  "type": "population_parameters_change",
  "time": 1000,
  "population": "pop1",
  "size": {"min": 5000, "max": 15000},  // Size range
  "growth_rate": {"min": 0.001, "max": 0.01}  // Growth rate range
}
```

### Migration Events

```json
{
  "type": "add_migration_rate_change",
  "time": {"min": 500, "max": 1500},  // Time range
  "source": "pop1",
  "dest": "pop2",
  "rate": {"min": 1e-5, "max": 1e-3}  // Rate range
}
```

### Population Split

```json
{
  "type": "split",
  "time": 2000,
  "source": "ancestral",
  "derived": ["pop1", "pop2"],
  "proportions": [0.5, 0.5]
}
```

### Mass Migration

```json
{
  "type": "mass_migration",
  "time": 1500,
  "source": "pop1",
  "dest": "pop2",
  "proportion": {"min": 0.1, "max": 0.3}  // Proportion range
}
```

---

## 5. Complete Examples

### Example 1: Simple Fixed Parameters

```json
{
  "name": "Fixed_Baseline",
  "description": "Baseline with all fixed parameters",
  
  "simulation_params": {
    "fasta_gz": "small_test_genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {
      "SNP": 50, "DEL": 20, "INS": 15, "INV": 10, "DUP": 5
    },
    "seed": "baseline"
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

### Example 2: Parameter Uncertainty

```json
{
  "name": "Uncertainty_Model",
  "description": "Model with parameter uncertainty",
  
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 2,
    "sv_distribution": {
      "SNP": {"min": 40, "max": 60},
      "DEL": {"min": 15, "max": 25},
      "INS": {"min": 10, "max": 20},
      "INV": {"min": 5, "max": 15},
      "DUP": {"min": 0, "max": 10}
    }
  },
  
  "evolutionary_params": {
    "mutation_rate": {"min": 5e-8, "max": 2e-7},
    "recombination_rate": {"min": 1e-9, "max": 1e-7},
    "generation_time": 25
  },
  
  "populations": [
    {"id": "pop1", "initial_size": {"min": 3000, "max": 8000}}
  ],
  
  "samples": [
    {"population": "pop1", "sample_size": 20}
  ]
}
```

### Example 3: Island Model with Migration

```json
{
  "name": "Island_Model",
  "description": "Five islands with migration",
  
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {
      "SNP": 40, "DEL": 20, "INS": 20, "INV": 15, "DUP": 5
    }
  },
  
  "evolutionary_params": {
    "mutation_rate": {"min": 1e-8, "max": 1e-6},
    "recombination_rate": 1e-8,
    "generation_time": 25
  },
  
  "populations": [
    {"id": "island1", "initial_size": 1000},
    {"id": "island2", "initial_size": {"min": 800, "max": 1200}},
    {"id": "island3", "initial_size": {"min": 800, "max": 1200}},
    {"id": "island4", "initial_size": 1000},
    {"id": "island5", "initial_size": {"min": 1500, "max": 2000}}
  ],
  
  "samples": [
    {"population": "island1", "sample_size": 5},
    {"population": "island2", "sample_size": 5},
    {"population": "island3", "sample_size": 5},
    {"population": "island4", "sample_size": 5},
    {"population": "island5", "sample_size": 5}
  ],
  
  "demographic_events": [
    {
      "type": "add_migration_rate_change",
      "time": 1000,
      "source": "island1",
      "dest": "island2",
      "rate": {"min": 5e-5, "max": 2e-4}
    }
  ]
}
```

---

## 6. Sampling Behavior

### Log-Uniform Sampling
Very small parameters (mutation/recombination rates < 0.01) use log-uniform sampling for better coverage across orders of magnitude.

### Uniform Sampling
Most parameters use standard uniform sampling between min and max.

### SV Distribution Normalization
After sampling, SV percentages are automatically normalized to sum to 100%.

---

## ✅ 7. Validation Checklist

Before using your model:

- ✔ All populations have `id` and `initial_size`
- ✔ Every sample references an existing population
- ✔ `simulation_params` includes `fasta_gz`, `chr_n`, and `sv_distribution`
- ✔ `evolutionary_params` includes `mutation_rate` and `recombination_rate`
- ✔ SV distributions sum to approximately 100% (will be normalized)
- ✔ All ranges have `min` < `max`
- ✔ Event times are non-negative

