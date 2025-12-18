#!/usr/bin/env python3
"""
Sample from parameter ranges in demographic models for MSpangepop

This script processes demographic models that contain parameter ranges and creates
individual configurations by randomly sampling from these ranges for each replicate.

Range notation in JSON:
- Fixed value: "mutation_rate": 1e-7
- Range: "mutation_rate": {"min": 1e-8, "max": 1e-6}

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
from typing import Any, Dict, Union

from io_handler import MSerror, MSsuccess, MScompute, MSwarning


# ============================================================================
# RANGE DETECTION AND SAMPLING
# ============================================================================

def is_range(value: Any) -> bool:
    """Check if a value is a range specification (dict with 'min' and 'max')."""
    return isinstance(value, dict) and 'min' in value and 'max' in value


def sample_from_range(range_spec: Dict[str, Any], param_name: str = "") -> float:
    """Sample a random value from a range specification."""
    import math
    
    min_val = range_spec['min']
    max_val = range_spec['max']
    
    if min_val > max_val:
        raise ValueError(f"Invalid range for {param_name}: min ({min_val}) > max ({max_val})")
    
    distribution = range_spec.get('distribution', 'auto').lower()
    
    if distribution == 'auto':
        if param_name in ['mutation_rate', 'recombination_rate'] and max_val < 0.01:
            distribution = 'log_uniform'
        else:
            distribution = 'uniform'
    
    if distribution in ['uniform', 'unif']:
        return random.uniform(min_val, max_val)
    
    elif distribution in ['log_uniform', 'loguniform', 'log']:
        if min_val <= 0:
            raise ValueError(f"Log-uniform requires positive values for {param_name}")
        log_min = math.log10(min_val)
        log_max = math.log10(max_val)
        return 10 ** random.uniform(log_min, log_max)
    
    elif distribution in ['normal', 'gaussian', 'norm']:
        mean = (min_val + max_val) / 2
        std = (max_val - min_val) / 6
        for _ in range(1000):
            sample = np.random.normal(mean, std)
            if min_val <= sample <= max_val:
                return float(sample)
        return random.uniform(min_val, max_val)
    
    else:
        MSwarning(f"Unknown distribution '{distribution}' for {param_name}, using uniform")
        return random.uniform(min_val, max_val)


def sample_dict_values(d: Dict, prefix: str = "") -> Dict:
    """Recursively sample from any range values in a dictionary."""
    result = {}
    for key, value in d.items():
        param_name = f"{prefix}_{key}" if prefix else key
        if is_range(value):
            result[key] = sample_from_range(value, param_name)
        elif isinstance(value, dict) and not is_range(value):
            result[key] = sample_dict_values(value, param_name)
        else:
            result[key] = value
    return result


def sample_sv_distribution(sv_ranges: Dict[str, Union[float, Dict]]) -> Dict[str, float]:
    """Sample from SV distribution ranges and normalize to sum to 100."""
    sampled = {}
    for sv_type, value in sv_ranges.items():
        if is_range(value):
            sampled[sv_type] = sample_from_range(value, f"sv_{sv_type}")
        else:
            sampled[sv_type] = value
    
    # Normalize to 100
    total = sum(sampled.values())
    if total > 0 and abs(total - 100) > 0.01:
        for sv_type in sampled:
            sampled[sv_type] = round((sampled[sv_type] / total) * 100, 2)
    
    return sampled


def sample_max_sv(max_sv_ranges: Dict[str, Union[int, Dict]]) -> Dict[str, int]:
    """Sample from max_sv ranges."""
    if max_sv_ranges is None:
        return None
    
    sampled = {}
    for sv_type, value in max_sv_ranges.items():
        if is_range(value):
            sampled[sv_type] = int(round(sample_from_range(value, f"max_sv_{sv_type}")))
        elif value is None:
            sampled[sv_type] = None
        else:
            sampled[sv_type] = int(value)
    
    return sampled


def process_demographic_ranges(demo_data: Dict, sample_name: str) -> Dict:
    """Process a demographic model and sample from any ranges."""
    sampled_demo = copy.deepcopy(demo_data)
    
    # Process evolutionary parameters
    if 'evolutionary_params' in sampled_demo:
        sampled_demo['evolutionary_params'] = sample_dict_values(
            sampled_demo['evolutionary_params'], 'evol'
        )
    
    # Process population parameters
    if 'populations' in sampled_demo:
        for pop in sampled_demo['populations']:
            if 'initial_size' in pop and is_range(pop['initial_size']):
                pop['initial_size'] = int(sample_from_range(pop['initial_size'], f"pop_{pop['id']}_size"))
    
    # Process demographic events
    if 'demographic_events' in sampled_demo:
        for i, event in enumerate(sampled_demo['demographic_events']):
            for key in ['time', 'size', 'growth_rate', 'proportion', 'rate']:
                if key in event and is_range(event[key]):
                    event[key] = sample_from_range(event[key], f"event_{i}_{key}")
    
    # Process migration matrix
    if 'migration_matrix' in sampled_demo:
        for migration in sampled_demo['migration_matrix']:
            if 'rate' in migration and is_range(migration['rate']):
                migration['rate'] = sample_from_range(migration['rate'], "migration_rate")
    
    # Process simulation_params
    if 'simulation_params' in sampled_demo:
        sim_params = sampled_demo['simulation_params']
        
        # SV distribution
        if 'sv_distribution' in sim_params:
            sim_params['sv_distribution'] = sample_sv_distribution(sim_params['sv_distribution'])
        
        # Max SV
        if 'max_sv' in sim_params:
            sim_params['max_sv'] = sample_max_sv(sim_params['max_sv'])
        
        # Other simple range params
        for key in ['minimal_sv_length']:
            if key in sim_params and is_range(sim_params[key]):
                sim_params[key] = int(sample_from_range(sim_params[key], key))
    
    return sampled_demo


# ============================================================================
# CONFIGURATION PROCESSING
# ============================================================================

def load_config(config_path: str) -> Dict:
    """Load the master configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_demographic_template(model_path: str) -> Dict:
    """Load a demographic model JSON file."""
    with open(model_path, 'r') as f:
        return json.load(f)


def expand_simulations(config: Dict) -> Dict:
    """Process samples and expand replicates with range sampling."""
    if 'samples' not in config:
        raise MSerror("No samples found in config")
    
    new_samples = {}
    expanded_dir = Path(".config/expanded_demographics")
    expanded_dir.mkdir(parents=True, exist_ok=True)
    
    for sample_base_name, sample_config in config['samples'].items():
        model_path = sample_config['model']
        replicates = sample_config.get('replicates', 1)
        
        MScompute(f"Model: {model_path}, Replicates: {replicates}")
        
        demo_template = load_demographic_template(model_path)
        
        for rep in range(1, replicates + 1):
            sample_name = f"{sample_base_name}_rep{rep}" if replicates > 1 else sample_base_name
            
            # Sample from all ranges
            sampled_demo = process_demographic_ranges(demo_template, sample_name)
            sampled_demo['name'] = sample_name
            
            # Save expanded demographic file
            demo_output_path = expanded_dir / f"{sample_name}_demographic.json"
            with open(demo_output_path, 'w') as f:
                json.dump(sampled_demo, f, indent=2)
            
            # Minimal YAML config - just point to demographic file
            sim_params = sampled_demo.get('simulation_params', {})
            new_samples[sample_name] = {
                'demographic_file': str(demo_output_path),
                'fasta_gz': sim_params.get('fasta_gz', 'reference.fa.gz'),
                'chr_n': sim_params.get('chr_n', 1),
            }
    
    MSsuccess(f"Generated {len(new_samples)} sample configurations")
    return new_samples


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    config_path = sys.argv[1] if len(sys.argv) > 1 else ".config/masterconfig.yaml"
    
    if len(sys.argv) > 2:
        random.seed(int(sys.argv[2]))
    
    MScompute(f"Loading config from: {config_path}")
    
    if not os.path.exists(config_path):
        raise MSerror(f"Config file not found: {config_path}")
    
    config = load_config(config_path)
    expanded_samples = expand_simulations(config)

    expanded_config = copy.deepcopy(config)
    expanded_config["samples"] = expanded_samples

    output_path = ".config/expanded_config.yaml"
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