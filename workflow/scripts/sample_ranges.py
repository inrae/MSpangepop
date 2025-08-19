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

def sample_from_range(range_spec: Dict[str, float], param_name: str = "") -> float:
    """
    Sample a random value from a range specification.
    
    Args:
        range_spec: Dictionary with 'min' and 'max' keys
        param_name: Name of the parameter (for logging)
    
    Returns:
        Random value between min and max
    """
    min_val = range_spec['min']
    max_val = range_spec['max']
    
    if min_val > max_val:
        raise ValueError(f"Invalid range for {param_name}: min ({min_val}) > max ({max_val})")
    
    # Use log-uniform distribution for very small values (like mutation rates)
    if param_name in ['mutation_rate', 'recombination_rate'] and max_val < 0.01:
        # Log-uniform sampling for rates
        import math
        log_min = math.log10(min_val)
        log_max = math.log10(max_val)
        log_sample = random.uniform(log_min, log_max)
        return 10 ** log_sample
    else:
        # Uniform sampling for other parameters
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

    
    # Process population parameters
    if 'populations' in sampled_demo:
        for pop in sampled_demo['populations']:
            if 'initial_size' in pop and is_range(pop['initial_size']):
                sampled_size = int(sample_from_range(pop['initial_size'], f"pop_{pop['id']}_size"))
                pop['initial_size'] = sampled_size
    
    # Process demographic events
    if 'demographic_events' in sampled_demo:
        for event in sampled_demo['demographic_events']:
            for key in ['time', 'size', 'growth_rate', 'proportion', 'rate']:
                if key in event and is_range(event[key]):
                    sampled_value = sample_from_range(event[key], f"event_{key}")
                    event[key] = sampled_value
    
    # Process migration matrix
    if 'migration_matrix' in sampled_demo:
        for migration in sampled_demo['migration_matrix']:
            if 'rate' in migration and is_range(migration['rate']):
                sampled_rate = sample_from_range(migration['rate'], "migration_rate")
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
    
    # Write output
    if expanded_samples:
        output_path = write_expanded_config(config, expanded_samples)
    else:
        output_path = config_path
    
    return output_path

if __name__ == "__main__":
    main()