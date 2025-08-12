#!/usr/bin/env python3
"""
Expand parameter sweeps into individual sample configurations for MSpangepop

This script takes a master configuration file containing parameter sweep definitions
and expands them into individual sample configurations. Each combination of parameters,
SV distributions, and replicates becomes a separate sample configuration.

MAIN WORKFLOW:
1. Load master configuration file (YAML)
2. For each sweep_sample in the config:
   - Load base demographic model (JSON)
   - Generate all combinations of demographic parameters
   - Generate all combinations with SV distributions  
   - Create replicates if specified
   - Save individual demographic files for each combination
3. Write expanded configuration with all individual samples

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

# ============================================================================
# CONFIGURATION LOADING AND VALIDATION
# ============================================================================

def load_config(config_path):
    """
    Load the master configuration file containing sweep definitions.
    """
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def validate_sv_distribution(sv_dist):
    """
    Validate that structural variant (SV) distribution percentages sum to 100.
    Also ensures all required SV types are present (sets missing ones to 0).
    
    Args:
        sv_dist (dict): Dictionary with SV types as keys and percentages as values
                       e.g., {'SNP': 70, 'DEL': 20, 'INS': 10}

    """
    total = sum(sv_dist.values())
    if abs(total - 100) > 0.01:
        raise ValueError(f"SV distribution sums to {total}, must be 100")
    
    # Ensure all required SV types are present (set missing to 0)
    required_types = {'SNP', 'DEL', 'INS', 'INV', 'DUP'}
    for sv_type in required_types:
        if sv_type not in sv_dist:
            sv_dist[sv_type] = 0
    
    return True

# ============================================================================
# STRING FORMATTING FOR FILENAMES
# ============================================================================

def format_param_value(value):
    """
    Format parameter values for use in filenames.
    Converts scientific notation to a filename-safe format.
    
    Args:
        value: Parameter value (float, int, or string)
        
    Returns:
        str: Filename-safe representation of the value
        
    Example:
        1e-4 -> "1e-4" -> "1em4" (replacing - with m for minus)
        0.001 -> "1e-3" -> "1em3"
    """
    if isinstance(value, float):
        return f"{value:.0e}".replace('+', '').replace('-', 'm')
    return str(value)

def format_sv_distribution(sv_dist):
    """
    Create a compact string representation of SV distribution for filenames.
    Uses single-letter abbreviations for each SV type.
    
    Args:
        sv_dist (dict): SV distribution dictionary
        
    Returns:
        str: Compact representation like "S70D20I10" for SNP:70%, DEL:20%, INS:10%
    """
    parts = []
    abbreviations = {'SNP': 'S', 'DEL': 'D', 'INS': 'I', 'INV': 'V', 'DUP': 'U'}
    
    for sv_type in ['SNP', 'DEL', 'INS', 'INV', 'DUP']:
        if sv_type in sv_dist and sv_dist[sv_type] > 0:
            parts.append(f"{abbreviations[sv_type]}{sv_dist[sv_type]}")
    
    return "".join(parts)

# ============================================================================
# NUMERIC PARAMETER PROCESSING
# ============================================================================

def ensure_numeric(value, param_name):
    """
    Ensure that numeric parameters are stored as proper numbers, not strings.
    Handles scientific notation and converts string representations to floats.
    
    This is important because YAML/JSON loading can sometimes read numeric
    values as strings, which can cause issues in downstream processing.
    """
    # Parameters that should always be numeric
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
    
    # Check if this parameter should be numeric based on its name
    if any(numeric_key in param_name.lower() for numeric_key in numeric_params):
        if isinstance(value, str):
            try:
                return float(value)  # Convert string to float
            except ValueError:
                MSwarning(f"Warning: Could not convert {param_name}={value} to number")
                return value
        elif isinstance(value, (int, float)):
            # Keep as float if it has decimal places or scientific notation
            return float(value) if '.' in str(value) or 'e' in str(value).lower() else value
    
    return value

def process_demographic_params(demo_copy, param_names, param_values):
    """
    Update demographic parameters in a demographic model, ensuring proper numeric types.
    Handles both top-level parameters and nested parameters within demographic events.
    """
    # Update top-level swept parameters
    for param_name, param_value in zip(param_names, param_values):
        if param_name in demo_copy:
            demo_copy[param_name] = ensure_numeric(param_value, param_name)
    
    # Ensure all top-level numeric fields are properly typed
    for key in ['mutation_rate', 'recombination_rate', 'generation_time']:
        if key in demo_copy:
            demo_copy[key] = ensure_numeric(demo_copy[key], key)
    
    # Process population parameters
    if 'populations' in demo_copy:
        for pop in demo_copy['populations']:
            if 'initial_size' in pop:
                pop['initial_size'] = ensure_numeric(pop['initial_size'], 'initial_size')
    
    # Process demographic events (size changes, splits, mergers, etc.)
    if 'demographic_events' in demo_copy:
        for event in demo_copy['demographic_events']:
            for key in ['time', 'size', 'growth_rate', 'proportion', 'rate']:
                if key in event:
                    event[key] = ensure_numeric(event[key], key)
    
    # Process migration matrix parameters
    if 'migration_matrix' in demo_copy:
        for migration in demo_copy['migration_matrix']:
            for key in ['time', 'rate']:
                if key in migration:
                    migration[key] = ensure_numeric(migration[key], key)
    
    return demo_copy

# ============================================================================
# MAIN SWEEP EXPANSION LOGIC
# ============================================================================

def expand_sweep_samples(config):
    """
    Main function to expand parameter sweeps into individual sample configurations.
    
    This is the core function that processes each sweep_sample definition and creates
    all possible combinations of:
    - Demographic parameters (if specified)
    - SV distributions
    - Replicates (if specified)
    
    For each combination, it:
    1. Creates a modified demographic JSON file
    2. Generates a sample configuration entry
    3. Assigns appropriate names and file paths
    """
    new_samples = {}
    
    # Create directory for expanded demographic files
    sweep_dir = Path(".config/expanded_demographics")
    sweep_dir.mkdir(parents=True, exist_ok=True)
    
    if 'sweep_samples' not in config:
        return new_samples
    
    MScompute("EXPANDING PARAMETER SWEEPS")
    
    # Process each sweep definition
    for sweep_name, sweep_config in config['sweep_samples'].items():
        MScompute(f"Processing sweep: {sweep_name}")
        
        # ===== STEP 1: Load base demographic model =====
        base_demo_path = sweep_config['base_demographic_file']
        MScompute(f"   Base model: {base_demo_path}")
        
        with open(base_demo_path, 'r') as f:
            base_demo = json.load(f)
        
        # ===== STEP 2: Process demographic parameter sweeps =====
        demo_sweep_params = sweep_config.get('sweep_params', {})
        
        # Convert all parameter values to proper numeric types
        for param_name, values in demo_sweep_params.items():
            demo_sweep_params[param_name] = [ensure_numeric(v, param_name) for v in values]
        
        # ===== STEP 3: Get SV distributions and other settings =====
        sv_distributions = sweep_config.get('sv_sweep', [])
        replicates = sweep_config.get('replicates', 1)
        seed_value = sweep_config.get('seed', None)  # Optional seed for reproducibility
        
        # Validate that SV distributions are provided
        if not sv_distributions:
            raise MSerror("Please specify a sv_distribution")
        
        # Validate all SV distributions sum to 100%
        for sv_dist in sv_distributions:
            try:
                validate_sv_distribution(sv_dist)
            except ValueError as e:
                MSwarning(f"{e} in distribution {sv_dist}")
                continue
        
        # ===== STEP 4: Generate all demographic parameter combinations =====
        demo_combinations = []
        if demo_sweep_params:
            param_names = list(demo_sweep_params.keys())
            param_values = [demo_sweep_params[k] for k in param_names]
            
            # Create all possible combinations using itertools.product
            for combination in itertools.product(*param_values):
                # Create filename-safe string representation of this combination
                param_str = "_".join([f"{k.replace('_rate', '')}{format_param_value(v)}" 
                                     for k, v in zip(param_names, combination)])
                demo_combinations.append((param_names, combination, param_str))
        else:
            # No demographic parameters to sweep
            demo_combinations.append(([], [], ""))
        
        # Calculate total number of configurations that will be generated
        total_runs = len(demo_combinations) * len(sv_distributions) * replicates
        MScompute(f"   Will generate {total_runs} configurations")
        
        # ===== STEP 5: Generate all combinations =====
        run_count = 0
        
        # Loop through demographic parameter combinations
        for demo_params, demo_values, demo_str in demo_combinations:
            
            # Loop through SV distributions
            for sv_idx, sv_distribution in enumerate(sv_distributions):
                sv_str = format_sv_distribution(sv_distribution)
                
                # Combine parameter strings for filename
                if demo_str and sv_str:
                    combined_str = f"{demo_str}_{sv_str}"
                elif demo_str:
                    combined_str = demo_str
                elif sv_str:
                    combined_str = sv_str
                else:
                    combined_str = "baseline"
                
                # Generate replicates
                for rep in range(1, replicates + 1):
                    run_count += 1
                    
                    # Create unique sample name
                    if replicates > 1:
                        sample_name = f"{sweep_name}_{combined_str}_rep{rep}"
                    else:
                        sample_name = f"{sweep_name}_{combined_str}"
                    
                    # ===== STEP 6: Create modified demographic file =====
                    demo_copy = copy.deepcopy(base_demo)
                    
                    # Update demographic parameters with proper numeric typing
                    demo_copy = process_demographic_params(demo_copy, demo_params, demo_values)
                    
                    # Update the model name to include sample name
                    demo_copy['name'] = f"{demo_copy.get('name', 'Model')}_{sample_name}"
                    
                    # Save expanded demographic file
                    demo_path = sweep_dir / f"{sample_name}_demographic.json"
                    with open(demo_path, 'w') as f:
                        json.dump(demo_copy, f, indent=2)
                    
                    # ===== STEP 7: Create sample configuration entry =====
                    new_samples[sample_name] = {
                        'fasta_gz': sweep_config['fasta_gz'],
                        'chr_n': sweep_config['chr_n'], 
                        'demographic_file': str(demo_path),
                        'sv_distribution': sv_distribution,
                    }
                    
                    # Add seed if specified (same seed for all expanded configs)
                    if seed_value is not None:
                        new_samples[sample_name]['seed'] = seed_value
                    
                    # Copy any additional parameters from sweep config
                    for key in ['model', 'minimal_sv_length', 'readable_json']:
                        if key in sweep_config:
                            new_samples[sample_name][key] = sweep_config[key]
                    
                    # Log progress for this configuration
                    sv_summary = ", ".join([f"{k}:{v}%" for k, v in sv_distribution.items() if v > 0])
                    MScompute(f"   [{run_count}/{total_runs}] {sample_name}: {sv_summary}")

    MSsuccess(f"Expanded {len(new_samples)} sample configurations")
    return new_samples

# ============================================================================
# OUTPUT WRITING
# ============================================================================

def write_expanded_config(config, expanded_samples, output_path=".config/expanded_config.yaml"):
    """
    Write the expanded configuration to a new YAML file.
    Merges the original config with the newly expanded samples.
    """
    # Create a deep copy to avoid modifying the original config
    expanded_config = copy.deepcopy(config)
    
    # Initialize samples section if it doesn't exist
    if 'samples' not in expanded_config:
        expanded_config['samples'] = {}
    
    # Add all expanded samples to the configuration
    expanded_config['samples'].update(expanded_samples)
    
    # Write to file with readable formatting
    with open(output_path, 'w') as f:
        yaml.dump(expanded_config, f, default_flow_style=False, sort_keys=False)
    
    MSsuccess(f"Expanded configuration written to: {output_path}")
    return output_path

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    # Parse command line arguments
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = ".config/masterconfig.yaml"
    
    MScompute(f"Loading config from: {config_path}")
    
    # Check that config file exists
    if not os.path.exists(config_path):
        MSerror(f"Config file not found: {config_path}")
        sys.exit(1)
    
    # Load and process configuration
    config = load_config(config_path)
    expanded_samples = expand_sweep_samples(config)
    
    # Write output
    if expanded_samples:
        output_path = write_expanded_config(config, expanded_samples)
        MSsuccess(f"Successfully expanded {len(expanded_samples)} samples")
    else:
        output_path = config_path
        MSwarning("No sweep samples found in configuration")
    
    return output_path

if __name__ == "__main__":
    main()