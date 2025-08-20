#  Creating a Demographic Model JSON for MSpangepop

## 📝 1. File Structure Overview

A complete demographic model now contains:
- **Simulation parameters** (genome file, chromosomes, SV distribution)
- **Evolutionary parameters** (mutation rate, recombination rate)
- **Population structure** (populations, samples, events)

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
    "recombination_rate": {"min": 1e-9, "max": 1e-7},  // Range with auto distribution
    "generation_time": 25                         // Fixed value
  }
}
```

### Range Syntax with Distribution Control

Use `{"min": value, "max": value, "distribution": "type"}` to specify both range and distribution:

```json
{
  "evolutionary_params": {
    "mutation_rate": {
      "min": 1e-8,
      "max": 1e-6,
      "distribution": "log_uniform"  // Explicitly use log-uniform
    },
    "recombination_rate": {
      "min": 1e-9,
      "max": 1e-7,
      "distribution": "normal"  // Sample from truncated normal
    }
  },
  
  "populations": [
    {
      "id": "pop1",
      "initial_size": {
        "min": 1000,
        "max": 10000,
        "distribution": "truncated_normal",  // Exponential distribution
        "mean": 2000,  // Optional parameter
        "std": 200
      }
    }
  ]
}
```

### Available Distributions

| Distribution | Keywords | Description | Extra Parameters |
|-------------|----------|-------------|------------------|
| **Uniform** | `"uniform"`, `"unif"` | Equal probability across range | None |
| **Log-uniform** | `"log_uniform"`, `"loguniform"`, `"log"` | Uniform in log space (good for rates) | None |
| **Normal** | `"normal"`, `"gaussian"`, `"norm"` | Bell curve, truncated to range | None |
| **Truncated Normal** | `"truncated_normal"`, `"truncnorm"` | Normal with custom mean/std | `mean`, `std` |
| **Auto** | `"auto"` or omitted | Auto-selects based on parameter type | None |

### Distribution Examples

#### Log-uniform for mutation rates (spans orders of magnitude)
```json
"mutation_rate": {
  "min": 1e-9,
  "max": 1e-6,
  "distribution": "log_uniform"
}
```

#### Normal distribution for population sizes
```json
"initial_size": {
  "min": 1000,
  "max": 10000,
  "distribution": "normal"  // Mean=(min+max)/2, truncated to range
}
```

#### Custom truncated normal
```json
"initial_size": {
  "min": 1000,
  "max": 10000,
  "distribution": "truncated_normal",
  "mean": 3000,  // Custom mean (closer to min)
  "std": 1500    // Custom standard deviation
}
```

### Auto Distribution Selection

When `"distribution"` is omitted or set to `"auto"`:
- **Mutation/recombination rates** (< 0.01): Uses `log_uniform`
- **All other parameters**: Uses `uniform`

### SV Distribution Ranges

Each SV type can have a range with distribution. The values are normalized to 100% after sampling:

```json
"sv_distribution": {
  "SNP": {
    "min": 30,
    "max": 50,
    "distribution": "normal"  // More likely to be near 40
  },
  "DEL": {"min": 15, "max": 25},  // Auto = uniform
  "INS": {"min": 10, "max": 20},
  "INV": {"min": 10, "max": 20},
  "DUP": {"min": 5, "max": 15}
}
```

---

## 📦 3. Parameter Sections

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

## 📋 5. Complete Examples

### Example 1: Simple Model with Distribution Control

```json
{
  "name": "Simple_Distribution_Example",
  "description": "Basic model showing distribution specifications",
  
  "simulation_params": {
    "fasta_gz": "small_test_genome.fa.gz",
    "chr_n": 1,
    "sv_distribution": {
      "SNP": 50, "DEL": 20, "INS": 15, "INV": 10, "DUP": 5
    },
    "seed": "baseline"
  },
  
  "evolutionary_params": {
    "mutation_rate": {
      "min": 1e-8,
      "max": 1e-6,
      "distribution": "log_uniform"  // Good for rates
    },
    "recombination_rate": 1e-8,  // Fixed
    "generation_time": 25
  },
  
  "populations": [
    {
      "id": "pop1",
      "initial_size": {
        "min": 3000,
        "max": 7000,
        "distribution": "normal"  // Bell curve around 5000
      }
    }
  ],
  
  "samples": [
    {"population": "pop1", "sample_size": 10}
  ]
}
```

### Example 2: Parameter Uncertainty with Custom Distributions

```json
{
  "name": "Advanced_Uncertainty_Model",
  "description": "Model demonstrating different distribution types",
  
  "simulation_params": {
    "fasta_gz": "genome.fa.gz",
    "chr_n": 2,
    "sv_distribution": {
      "SNP": {
        "min": 40,
        "max": 60,
        "distribution": "normal"  // Peak around 50%
      },
      "DEL": {"min": 15, "max": 25},  // Uniform (default)
      "INS": {"min": 10, "max": 20},
      "INV": {
        "min": 5,
        "max": 15,
        "distribution": "normal"
      },
      "DUP": {"min": 0, "max": 10}
    }
  },
  
  "evolutionary_params": {
    "mutation_rate": {
      "min": 1e-9,
      "max": 1e-6,
      "distribution": "log_uniform"  
    },
    "recombination_rate": {
      "min": 1e-9,
      "max": 1e-7,
      "distribution": "truncated_normal",
      "mean": 5e-8,  // Peak around this value
      "std": 2e-8
    },
    "generation_time": 25  // Fixed
  },
  
  "populations": [
    {
      "id": "pop1",
      "initial_size": {
        "min": 1000,
        "max": 10000
      }
    },
    {
      "id": "pop2",
      "initial_size": {
        "min": 500,
        "max": 5000,
        "distribution": "normal"
      }
    }
  ],
  
  "samples": [
    {"population": "pop1", "sample_size": 10},
    {"population": "pop2", "sample_size": 10}
  ],
  
  "demographic_events": [
    {
      "type": "mass_migration",
      "time": 1000,
      "source": "pop2",
      "dest": "pop1",
      "proportion": {
        "min": 0.1,
        "max": 0.5
      }
    }
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

### SV Distribution Normalization
After sampling all SV types, percentages are automatically normalized to sum to 100%:

1. Sample each type according to its distribution
2. Sum all sampled values
3. Scale each value by (100 / sum)

Example:
- SNP samples 45 (from normal distribution)
- DEL samples 22 (from uniform)
- INS samples 18 (from uniform)
- INV samples 20 (from uniform)
- DUP samples 10 (from uniform)
- Total: 115 → Each scaled by 100/115

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
