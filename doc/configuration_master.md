# Configuring the Pipeline - Master Configuration

## Overview

The master configuration file (`.config/masterconfig.yaml`) is the entry point for running MSpangepop simulations. It defines which demographic models to run and how many replicates to generate for each.

## Configuration Structure

The configuration file uses a simple YAML format with two main sections:

```yaml
# Named sample configurations
samples:
  sample_name:
    model: "path/to/demographic_model.json"
    replicates: int

# Global settings
output_dir: "results/"
memory_multiplier: 1
```

## Samples Section

Each entry in the `samples` section defines a simulation run:

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `model` | string | Path to the demographic model JSON file | `"simulation_data/Panmictic_Model.json"` |
| `replicates` | integer | Number of simulation replicates to run | `10` |

### Sample Naming Convention

- **Single replicate** (`replicates: 1`): The sample name is used as-is
  - Config: `baseline` → Output directory: `results/baseline/`
  
- **Multiple replicates** (`replicates: > 1`): A `_repN` suffix is added
  - Config: `test_run` with `replicates: 3`
  - Output directories: `results/test_run_rep1/`, `results/test_run_rep2/`, `results/test_run_rep3/`

### Example Configurations

#### Simple Test Run
```yaml
samples:
  quick_test:
    model: "simulation_data/Test_Model.json"
    replicates: 1
```

#### Multiple Models
```yaml
samples:
  panmictic_baseline:
    model: "simulation_data/Panmictic_Model.json"
    replicates: 1
    
  island_analysis:
    model: "simulation_data/Island_Model.json"
    replicates: 5
    
  complex_demography:
    model: "simulation_data/Complex_Model.json"
    replicates: 10
```

#### Parameter Exploration with Ranges
```yaml
samples:
  # Fixed parameters - good for baseline
  fixed_params:
    model: "simulation_data/Fixed_Model.json"
    replicates: 1
    
  # Range parameters - explores uncertainty
  parameter_uncertainty:
    model: "simulation_data/Model_With_Ranges.json"
    replicates: 100  # Each replicate samples different values from ranges
```

## Global Settings

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `output_dir` | string | Base directory for all simulation outputs | `"results/"` |
| `memory_multiplier` | float | Memory scaling factor for cluster jobs (increase if OOM errors) | `1` |
| `succint` | bool | If True, skip visialization rule | `False` |


## Complete Example

```yaml
# .config/masterconfig.yaml

# Named sample configurations
samples:
  # Baseline simulation with fixed parameters
  baseline:
    model: "simulation_data/Fixed_Baseline.json"
    replicates: 1
  
  # Test different parameter values
  mutation_rate_test:
    model: "simulation_data/Mutation_Range_Model.json"
    replicates: 20  # Each rep samples different mutation rate
  
  # Compare different demographic scenarios
  panmictic:
    model: "simulation_data/Panmictic_Model.json"
    replicates: 10
    
  island:
    model: "simulation_data/Island_Model.json"
    replicates: 10
  
  # Production run with uncertainty quantification
  production_analysis:
    model: "simulation_data/Production_Model_Ranges.json"
    replicates: 100

# Global settings
output_dir: "results/"
memory_multiplier: 1.5  
succint: False
```

## How Replicates Work with Parameter Ranges

When your demographic model contains parameter ranges (see [Demographic Model Configuration](configuration_model.md)):

1. **Fixed parameters**: All replicates use identical values
2. **Range parameters**: Each replicate randomly samples from the specified ranges
3. **Mixed models**: Fixed parameters stay constant, ranges are sampled

Example: If `Model_With_Ranges.json` contains:
```json
{
  "evolutionary_params": {
    "mutation_rate": {"min": 1e-8, "max": 1e-6},  // Range
    "recombination_rate": 1e-8,                    // Fixed
    "generation_time": 25                          // Fixed
  }
}
```

With `replicates: 10`, you'll get:
- 10 simulations with the SAME recombination rate (1e-8) and generation time (25)
- 10 simulations with DIFFERENT mutation rates (randomly sampled between 1e-8 and 1e-6)

## Workflow Execution

When you run the workflow:

1. **Expansion phase**: The `sample_ranges.py` script processes your configuration
   - Creates individual configurations for each replicate
   - Samples values from any parameter ranges
   - Saves expanded demographic files to `.config/expanded_demographics/`

2. **Simulation phase**: Each expanded sample runs through the full pipeline
   - Coalescent simulation with msprime
   - Variant generation
   - Graph construction

3. **Output organization**: Results are saved to named directories
   - Single replicate: `results/sample_name/`
   - Multiple replicates: `results/sample_name_rep1/`, `results/sample_name_rep2/`, etc.


## Troubleshooting

**Q: Where are the simulation parameters defined?**  
A: All simulation parameters (mutation rate, population sizes, SV distributions, etc.) are  defined in the demographic model JSON files, not in masterconfig.yaml.

**Q: How do I run the same model with different parameters?**  
A: Use parameter ranges in your demographic model JSON and set `replicates` > 1.

