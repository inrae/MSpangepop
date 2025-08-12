#!/usr/bin/env python3
"""
Expand parameter sweeps into individual sample configurations for MSpangepop
Author: Lucien Piat
"""

import yaml
import json
import itertools
import os
import sys
from pathlib import Path
import copy
from io_handler import MSerror, MSsuccess, MScompute, MSwarning

def load_config(config_path):
    """Load the master config file"""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def format_param_value(value):
    """Format parameter value for filename"""
    if isinstance(value, float):
        return f"{value:.0e}".replace('+', '').replace('-', 'm')
    return str(value)

def format_sv_distribution(sv_dist):
    """Create a short string representation of SV distribution for filenames"""
    parts = []
    abbreviations = {'SNP': 'S', 'DEL': 'D', 'INS': 'I', 'INV': 'V', 'DUP': 'U'}
    for sv_type in ['SNP', 'DEL', 'INS', 'INV', 'DUP']:
        if sv_type in sv_dist and sv_dist[sv_type] > 0:
            parts.append(f"{abbreviations[sv_type]}{sv_dist[sv_type]}")
    return "".join(parts)

def validate_sv_distribution(sv_dist):
    """Validate that SV distribution sums to 100"""
    total = sum(sv_dist.values())
    if abs(total - 100) > 0.01:
        raise ValueError(f"SV distribution sums to {total}, must be 100")
    
    required_types = {'SNP', 'DEL', 'INS', 'INV', 'DUP'}
    for sv_type in required_types:
        if sv_type not in sv_dist:
            sv_dist[sv_type] = 0
    
    return True

def ensure_numeric(value, param_name):
    """
    Ensure that numeric parameters are stored as numbers, not strings.
    Handles scientific notation properly.
    """
    # List of parameters that should always be numeric
    numeric_params = [
        'mutation_rate', 
        'recombination_rate', 
        'generation_time',
        'initial_size',
        'growth_rate',
        'proportion',
        'rate',
        'time'
    ]
    
    # Check if this parameter should be numeric
    if any(numeric_key in param_name.lower() for numeric_key in numeric_params):
        if isinstance(value, str):
            try:
                # Try to convert to float
                return float(value)
            except ValueError:
                MSwarning(f"Warning: Could not convert {param_name}={value} to number")
                return value
        elif isinstance(value, (int, float)):
            return float(value) if '.' in str(value) or 'e' in str(value).lower() else value
    
    return value

def process_demographic_params(demo_copy, param_names, param_values):
    """
    Update demographic parameters ensuring proper numeric types.
    Also handles nested parameters in demographic events.
    """
    # Update top-level parameters
    for param_name, param_value in zip(param_names, param_values):
        if param_name in demo_copy:
            demo_copy[param_name] = ensure_numeric(param_value, param_name)
    
    # Ensure all top-level numeric fields are properly typed
    for key in ['mutation_rate', 'recombination_rate', 'generation_time']:
        if key in demo_copy:
            demo_copy[key] = ensure_numeric(demo_copy[key], key)
    
    # Process populations
    if 'populations' in demo_copy:
        for pop in demo_copy['populations']:
            if 'initial_size' in pop:
                pop['initial_size'] = ensure_numeric(pop['initial_size'], 'initial_size')
    
    # Process demographic events
    if 'demographic_events' in demo_copy:
        for event in demo_copy['demographic_events']:
            for key in ['time', 'size', 'growth_rate', 'proportion', 'rate']:
                if key in event:
                    event[key] = ensure_numeric(event[key], key)
    
    # Process migration matrix
    if 'migration_matrix' in demo_copy:
        for migration in demo_copy['migration_matrix']:
            for key in ['time', 'rate']:
                if key in migration:
                    migration[key] = ensure_numeric(migration[key], key)
    
    return demo_copy

def expand_parameter_sweeps(config):
    """
    Expand parameter sweeps into individual sample configurations
    """
    
    new_samples = {}
    sweep_dir = Path(".config/expanded_demographics")
    sweep_dir.mkdir(parents=True, exist_ok=True)
    
    if 'parameter_sweeps' not in config:
        return new_samples
    
    MScompute("EXPANDING PARAMETER SWEEPS")
    
    for sweep_name, sweep_config in config['parameter_sweeps'].items():

        # Load base demographic file
        base_demo_path = sweep_config['base_demographic_file']
        MScompute(f"   Base model: {base_demo_path}")
        
        with open(base_demo_path, 'r') as f:
            base_demo = json.load(f)
        
        # Get sweep parameters and ensure they're numeric
        demo_sweep_params = sweep_config.get('sweep_params', {})
        
        # Convert parameter values to proper numeric types
        for param_name, values in demo_sweep_params.items():
            demo_sweep_params[param_name] = [ensure_numeric(v, param_name) for v in values]
        
        sv_distributions = sweep_config.get('sv_sweep', [])
        replicates = sweep_config.get('replicates', 1)
        
        # Handle seed - NEW LOGIC
        seed_value = sweep_config.get('seed', None)
        
        # If no SV sweep specified, use default
        if not sv_distributions:
            default_sv = sweep_config.get('sv_distribution', 
                                         {'SNP': 50, 'DEL': 20, 'INS': 20, 'INV': 10, 'DUP': 0})
            sv_distributions = [default_sv]
        
        # Validate all SV distributions
        for sv_dist in sv_distributions:
            try:
                validate_sv_distribution(sv_dist)
            except ValueError as e:
                MSwarning(f"{e} in distribution {sv_dist}")
                continue
        
        # Handle demographic parameters
        demo_combinations = []
        if demo_sweep_params:
            param_names = list(demo_sweep_params.keys())
            param_values = [demo_sweep_params[k] for k in param_names]
            
            for combination in itertools.product(*param_values):
                param_str = "_".join([f"{k.replace('_rate', '')}{format_param_value(v)}" 
                                     for k, v in zip(param_names, combination)])
                demo_combinations.append((param_names, combination, param_str))
        else:
            demo_combinations.append(([], [], ""))
        
        # Calculate total runs
        total_runs = len(demo_combinations) * len(sv_distributions) * replicates

        
        # Generate all combinations
        run_count = 0
        for demo_params, demo_values, demo_str in demo_combinations:
            for sv_idx, sv_distribution in enumerate(sv_distributions):
                sv_str = format_sv_distribution(sv_distribution)
                
                # Combine parameter strings
                if demo_str and sv_str:
                    combined_str = f"{demo_str}_{sv_str}"
                elif demo_str:
                    combined_str = demo_str
                elif sv_str:
                    combined_str = sv_str
                else:
                    combined_str = "baseline"
                
                # Add replicates
                for rep in range(1, replicates + 1):
                    run_count += 1
                    
                    # Create unique sample name
                    if replicates > 1:
                        sample_name = f"{sweep_name}_{combined_str}_rep{rep}"
                    else:
                        sample_name = f"{sweep_name}_{combined_str}"
                    
                    # Create modified demographic file
                    demo_copy = copy.deepcopy(base_demo)
                    
                    # Process demographic parameters with proper typing
                    demo_copy = process_demographic_params(demo_copy, demo_params, demo_values)
                    
                    # Update the name in the demographic file
                    demo_copy['name'] = f"{demo_copy.get('name', 'Model')}_{sample_name}"
                    
                    # Save expanded demographic file
                    demo_path = sweep_dir / f"{sample_name}_demographic.json"
                    with open(demo_path, 'w') as f:
                        json.dump(demo_copy, f, indent=2)
                    
                    # Create sample configuration
                    new_samples[sample_name] = {
                        'fasta_gz': sweep_config['fasta_gz'],
                        'chr_n': sweep_config['chr_n'],
                        'demographic_file': str(demo_path),
                        'sv_distribution': sv_distribution,
                    }
                    
                    # FIXED SEED HANDLING:
                    # Only add seed if it was specified in the sweep config
                    # Use the SAME seed for all expanded configs
                    if seed_value is not None:
                        new_samples[sample_name]['seed'] = seed_value
                    
                    # Copy any additional parameters
                    for key in ['model', 'minimal_sv_length', 'readable_json']:
                        if key in sweep_config:
                            new_samples[sample_name][key] = sweep_config[key]
                    
                    
                    
                    # Show SV distribution summary
                    sv_summary = ", ".join([f"{k}:{v}%" for k, v in sv_distribution.items() if v > 0])

    return new_samples

# Rest of the script remains the same...
def write_expanded_config(config, expanded_samples, output_path=".config/expanded_config.yaml"):
    """Write the expanded configuration to a new file"""
    
    expanded_config = copy.deepcopy(config)
    
    if 'samples' not in expanded_config:
        expanded_config['samples'] = {}
    
    expanded_config['samples'].update(expanded_samples)
    
    with open(output_path, 'w') as f:
        yaml.dump(expanded_config, f, default_flow_style=False, sort_keys=False)
    
    
    return output_path

def main():
    """Main function to expand sweeps"""
    
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = ".config/masterconfig.yaml"
    
    MScompute(f"Loading config from: {config_path}")
    
    if not os.path.exists(config_path):
        MSerror(f"Config file not found: {config_path}")
        sys.exit(1)
    
    config = load_config(config_path)
    
    expanded_samples = expand_parameter_sweeps(config)
    
    if expanded_samples:
        output_path = write_expanded_config(config, expanded_samples)
    else:
        output_path = config_path
    
    return output_path

if __name__ == "__main__":
    main()