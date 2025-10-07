#!/usr/bin/env python3
"""
Sample from parameter ranges in demographic models for MSpangepop

This script processes demographic models that contain parameter ranges and creates
individual configurations by randomly sampling from these ranges for each replicate.

Range notation in JSON:
- Fixed value: "mutation_rate": 1e-7
- Range: "mutation_rate": {"min": 1e-8, "max": 1e-6}
- SV ranges: "SNP": {"min": 30, "max": 50}

Author: Lucien Piat
"""

import yaml
import json
import os
import sys
import copy
import random
import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Union

from io_handler import MSerror, MSsuccess, MScompute, MSwarning


# ============================================================================
# RANGE DETECTION AND SAMPLING
# ============================================================================

def is_range(value: Any) -> bool:
    """
    Check if a value is a range specification.
    Ranges are dictionaries with 'min' and 'max' keys.
    """
    return isinstance(value, dict) and 'min' in value and 'max' in value

def sample_from_range(range_spec: Dict[str, Any], param_name: str = "") -> float:
    """
    Sample a random value from a range specification with support for different distributions.
    
    Args:
        range_spec: Dictionary with 'min', 'max', and optional 'distribution' keys
        param_name: Name of the parameter (for logging)
    
    Returns:
        Random value sampled according to specified distribution
    """
    import math
    
    min_val = range_spec['min']
    max_val = range_spec['max']
    
    if min_val > max_val:
        raise ValueError(f"Invalid range for {param_name}: min ({min_val}) > max ({max_val})")
    
    # Get distribution type (default based on parameter type)
    distribution = range_spec.get('distribution', 'auto').lower()
    
    # Auto-detect best distribution if not specified
    if distribution == 'auto':
        if param_name in ['mutation_rate', 'recombination_rate'] and max_val < 0.01:
            distribution = 'log_uniform'
        else:
            distribution = 'uniform'
    
    # Sample based on distribution type
    if distribution in ['uniform', 'unif']:
        # Uniform distribution
        return random.uniform(min_val, max_val)
    
    elif distribution in ['log_uniform', 'loguniform', 'log']:
        # Log-uniform distribution (good for rates spanning orders of magnitude)
        if min_val <= 0:
            raise ValueError(f"Log-uniform requires positive values for {param_name}")
        log_min = math.log10(min_val)
        log_max = math.log10(max_val)
        log_sample = random.uniform(log_min, log_max)
        return 10 ** log_sample
    
    elif distribution in ['normal', 'gaussian', 'norm']:
        # Normal distribution (truncated to range)
        # Use min/max as roughly 3 sigma bounds
        mean = (min_val + max_val) / 2
        std = (max_val - min_val) / 6  # 99.7% within bounds
        
        # Sample until we get a value within bounds (truncated normal)
        for _ in range(1000):  # Max attempts to avoid infinite loop
            sample = np.random.normal(mean, std)
            if min_val <= sample <= max_val:
                return float(sample)
        # Fallback to uniform if we can't get a valid sample
        MSwarning(f"Could not sample from truncated normal for {param_name}, using uniform")
        return random.uniform(min_val, max_val)
    
    elif distribution in ['truncated_normal', 'truncnorm']:
        # Truncated normal with custom mean and std
        mean = range_spec.get('mean', (min_val + max_val) / 2)
        std = range_spec.get('std', (max_val - min_val) / 4)
        
        # Use scipy if available, otherwise fallback
        try:
            from scipy import stats
            a = (min_val - mean) / std
            b = (max_val - mean) / std
            return float(stats.truncnorm.rvs(a, b, loc=mean, scale=std))
        except ImportError:
            # Fallback without scipy
            for _ in range(1000):
                sample = np.random.normal(mean, std)
                if min_val <= sample <= max_val:
                    return float(sample)
            return random.uniform(min_val, max_val)
    
    else:
        MSwarning(f"Unknown distribution '{distribution}' for {param_name}, using uniform")
        return random.uniform(min_val, max_val)

def sample_sv_distribution(sv_ranges: Dict[str, Union[float, Dict]]) -> Dict[str, float]:
    """
    Sample from SV distribution ranges and normalize to sum to 100.
    
    Args:
        sv_ranges: Dictionary where values can be fixed numbers or range specs
    
    Returns:
        Dictionary with sampled SV percentages that sum to 100
    """
    sampled = {}
    
    # Sample each SV type
    for sv_type, value in sv_ranges.items():
        if is_range(value):
            sampled[sv_type] = sample_from_range(value, f"sv_{sv_type}")
        else:
            sampled[sv_type] = value
    
    # Normalize to sum to 100
    total = sum(sampled.values())
    if abs(total - 100) > 0.01:
        for sv_type in sampled:
            sampled[sv_type] = (sampled[sv_type] / total) * 100
    
    # Round to reasonable precision
    for sv_type in sampled:
        sampled[sv_type] = round(sampled[sv_type], 2)
    
    return sampled

def process_demographic_ranges(demo_data: Dict, sample_name: str) -> Dict:
    """
    Process a demographic model containing ranges and sample random values.
    
    Args:
        demo_data: Demographic model data (may contain ranges)
        sample_name: Name for this sample (for logging)
    
    Returns:
        New demographic model with all ranges replaced by sampled values
    """
    # Deep copy to avoid modifying original
    sampled_demo = copy.deepcopy(demo_data)
    
    # Process evolutionary parameters
    if 'evolutionary_params' in sampled_demo:
        params = sampled_demo['evolutionary_params']
        for param_name, value in params.items():
            if is_range(value):
                sampled_value = sample_from_range(value, param_name)
                params[param_name] = sampled_value
                dist = value.get('distribution', 'auto')
                # Determine actual distribution used if auto
                if dist == 'auto':
                    if param_name in ['mutation_rate', 'recombination_rate'] and value['max'] < 0.01:
                        dist = 'log_uniform'
                    else:
                        dist = 'uniform'
    
    # Process population parameters
    if 'populations' in sampled_demo:
        for pop in sampled_demo['populations']:
            if 'initial_size' in pop and is_range(pop['initial_size']):
                sampled_size = int(sample_from_range(pop['initial_size'], f"pop_{pop['id']}_size"))
                dist = pop['initial_size'].get('distribution', 'uniform')
                pop['initial_size'] = sampled_size
    
    # Process demographic events
    if 'demographic_events' in sampled_demo:
        for i, event in enumerate(sampled_demo['demographic_events']):
            for key in ['time', 'size', 'growth_rate', 'proportion', 'rate']:
                if key in event and is_range(event[key]):
                    sampled_value = sample_from_range(event[key], f"event_{i}_{key}")
                    dist = event[key].get('distribution', 'uniform')
                    event[key] = sampled_value
    
    # Process migration matrix
    if 'migration_matrix' in sampled_demo:
        for migration in sampled_demo['migration_matrix']:
            if 'rate' in migration and is_range(migration['rate']):
                sampled_rate = sample_from_range(migration['rate'], "migration_rate")
                dist = migration['rate'].get('distribution', 'uniform')
                migration['rate'] = sampled_rate
    
    return sampled_demo

# ============================================================================
# CONFIGURATION PROCESSING
# ============================================================================

def load_config(config_path: str) -> Dict:
    """Load the master configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_demographic_template(model_path: str) -> Dict:
    """Load a demographic model JSON file (may contain ranges)."""
    with open(model_path, 'r') as f:
        return json.load(f)

def expand_simulations(config: Dict) -> Dict:
    """
    Process samples and expand replicates with range sampling.
    
    Args:
        config: Master configuration dictionary
    
    Returns:
        Expanded configuration with individual samples
    """
    if 'samples' not in config:
        raise MSerror("No samples found in config")
    
    new_samples = {}
    expanded_dir = Path(".config/expanded_demographics")
    expanded_dir.mkdir(parents=True, exist_ok=True)
    
    for sample_base_name, sample_config in config['samples'].items():
        model_path = sample_config['model']
        replicates = sample_config.get('replicates', 1)
        
        MScompute(f"Model: {model_path}, Replicates: {replicates}")
        
        # Load the demographic template (may contain ranges)
        demo_template = load_demographic_template(model_path)
        
        # Extract simulation parameters from the demographic file
        sim_params = demo_template.get('simulation_params', {})
        
        for rep in range(1, replicates + 1):
            # Create unique sample name based on the base name from config
            if replicates > 1:
                sample_name = f"{sample_base_name}_rep{rep}"
            else:
                sample_name = sample_base_name
            
            # Sample from ranges in demographic parameters
            sampled_demo = process_demographic_ranges(demo_template, sample_name)
            
            # Ensure simulation_params exists in sampled_demo
            if 'simulation_params' not in sampled_demo:
                sampled_demo['simulation_params'] = {}
            
            # Copy simulation parameters from template
            sampled_demo['simulation_params'].update(sim_params)
            
            # Sample from SV distribution ranges if present
            if 'sv_distribution' in sim_params:
                sv_dist = sim_params['sv_distribution']
                sampled_sv = sample_sv_distribution(sv_dist)
                # Update the sampled demo with the sampled SV distribution
                sampled_demo['simulation_params']['sv_distribution'] = sampled_sv
                
                # Log sampled SV distribution
                sv_str = ", ".join([f"{k}:{v:.1f}%" for k, v in sampled_sv.items() if v > 0])
            
            # Update the model name to include sample name
            sampled_demo['name'] = sample_name
            
            # Save the sampled demographic file
            demo_output_path = expanded_dir / f"{sample_name}_demographic.json"
            with open(demo_output_path, 'w') as f:
                json.dump(sampled_demo, f, indent=2)
            
            # Create sample configuration
            new_samples[sample_name] = {
                'demographic_file': str(demo_output_path),
                'fasta_gz': sampled_demo['simulation_params'].get('fasta_gz', 'reference.fa.gz'),
                'chr_n': sampled_demo['simulation_params'].get('chr_n', 1),
                'sv_distribution': sampled_demo['simulation_params'].get('sv_distribution', {}),
                'seed': sampled_demo['simulation_params'].get('seed'),
            }
            
            # Copy optional parameters
            for key in ['model', 'minimal_sv_length', 'readable_json']:
                if key in sampled_demo['simulation_params'] and sampled_demo['simulation_params'][key] is not None:
                    new_samples[sample_name][key] = sampled_demo['simulation_params'][key]
    
    MSsuccess(f"Generated {len(new_samples)} sample configurations")
    return new_samples

def write_expanded_config(config: Dict, expanded_samples: Dict, output_path: str = ".config/expanded_config.yaml") -> str:
    """
    Write the expanded configuration to a new YAML file.
    
    Args:
        config: Original configuration
        expanded_samples: Dictionary of expanded sample configurations
        output_path: Path to write the expanded config
    
    Returns:
        Path to the written file
    """
    expanded_config = copy.deepcopy(config)
    
    # Replace the original samples with expanded ones
    expanded_config['samples'] = expanded_samples
    
    with open(output_path, 'w') as f:
        yaml.dump(expanded_config, f, default_flow_style=False, sort_keys=False)
    
    return output_path

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    # Parse command line arguments
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = ".config/masterconfig.yaml"
    
    # Set random seed for reproducibility if desired
    if len(sys.argv) > 2:
        random.seed(int(sys.argv[2]))
    
    MScompute(f"Loading config from: {config_path}")
    
    # Check that config file exists
    if not os.path.exists(config_path):
        raise MSerror(f"Config file not found: {config_path}")
    
    # Load and process configuration
    config = load_config(config_path)
    expanded_samples = expand_simulations(config)

    # Build the expanded config in memory
    expanded_config = copy.deepcopy(config)
    expanded_config["samples"] = expanded_samples

    output_path = ".config/expanded_config.yaml"

    # --- Only write if content changed ---
    new_yaml = yaml.dump(expanded_config, default_flow_style=False, sort_keys=False)
    if os.path.exists(output_path):
        with open(output_path, "r") as f:
            old_yaml = f.read()
    else:
        old_yaml = None

    if old_yaml == new_yaml:
        MSsuccess(f"{output_path} is already up-to-date. No rewrite needed.")
    else:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write(new_yaml)
        MSsuccess(f"Wrote updated expanded config to {output_path}.")

    return output_path


if __name__ == "__main__":
    main()
