#!/usr/bin/env python3

"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage: Create a recap for MSpangepop simulation
"""

import os
import yaml
import json
import argparse
import socket
import datetime
import platform
from io_handler import MSsuccess, MScompute

os.environ['MPLCONFIGDIR'] = './.config/matplotlib'

import pkg_resources


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a recap file containing simulation parameters."
    )
    parser.add_argument("config_file", help="Path to the configuration YAML file")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("current_run", help="Name of the current run/sample")
    return parser.parse_args()


def load_demographic_params(demographic_file):
    """Load simulation parameters from demographic JSON file."""
    with open(demographic_file, 'r') as f:
        return json.load(f)


def main():
    args = parse_args()
    recap_file = os.path.join(args.output_dir, f"{args.current_run}_global_recap.txt")
    os.makedirs(args.output_dir, exist_ok=True)

    # Load YAML config to get demographic file path
    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    sample_config = config.get("samples", {}).get(args.current_run, {})
    demographic_file = sample_config.get("demographic_file")
    
    # Load all params from demographic file
    demo_data = load_demographic_params(demographic_file) if demographic_file else {}
    sim_params = demo_data.get('simulation_params', {})
    evol_params = demo_data.get('evolutionary_params', {})

    # Global info
    machine_name = socket.gethostname()
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    python_version = platform.python_version()
    os_info = f"{platform.system()} {platform.release()}"
    current_dir = os.getcwd()

    installed_packages = {
        dist.project_name: dist.version
        for dist in pkg_resources.working_set
    }

    with open(recap_file, "w") as f:
        f.write(f"\U0001F539 MSpangepop Recap File for: {args.current_run}\n\n")

        f.write("\U0001F539 Global Run Information:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Machine Name       : {machine_name}\n")
        f.write(f"Date/Time          : {current_date}\n")
        f.write(f"Python Version     : {python_version}\n")
        f.write(f"Operating System   : {os_info}\n")
        f.write(f"Current Directory  : {current_dir}\n\n")

        f.write("\U0001F539 Sample Configuration:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Sample Name         : {args.current_run}\n")
        f.write(f"Demographic File    : {demographic_file}\n")
        f.write(f"FASTA File          : {sim_params.get('fasta_gz', 'N/A')}\n")
        f.write(f"Chromosome Count    : {sim_params.get('chr_n', 'N/A')}\n")
        f.write(f"Seed                : {sim_params.get('seed', 'None (random)')}\n")
        f.write(f"Minimal SV Length   : {sim_params.get('minimal_sv_length', 1)}\n")
        f.write(f"Readable JSON       : {sim_params.get('readable_json', False)}\n\n")

        # Evolutionary parameters
        if evol_params:
            f.write("\U0001F539 Evolutionary Parameters:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mutation Rate       : {evol_params.get('mutation_rate', 'N/A')}\n")
            f.write(f"Recombination Rate  : {evol_params.get('recombination_rate', 'N/A')}\n")
            f.write(f"Generation Time     : {evol_params.get('generation_time', 'N/A')}\n\n")

        # SV Distribution
        sv_dist = sim_params.get('sv_distribution', {})
        if sv_dist:
            f.write("\U0001F539 SV Distribution (%):\n")
            f.write("-" * 40 + "\n")
            for sv_type, percentage in sv_dist.items():
                f.write(f"  {sv_type:10}: {percentage}%\n")
            f.write("\n")

        # Max SV
        max_sv = sim_params.get('max_sv')
        if max_sv:
            f.write("\U0001F539 Max SV Counts (global limits):\n")
            f.write("-" * 40 + "\n")
            total_max = 0
            for sv_type, count in max_sv.items():
                if count is not None:
                    f.write(f"  {sv_type:10}: {count}\n")
                    total_max += count
                else:
                    f.write(f"  {sv_type:10}: unlimited\n")
            f.write(f"  {'TOTAL':10}: {total_max}\n\n")

        # SV Length Files
        sv_length_files = sim_params.get('sv_length_files')
        if sv_length_files:
            f.write("\U0001F539 Custom SV Length Distribution Files:\n")
            f.write("-" * 40 + "\n")
            for sv_type, file_path in sv_length_files.items():
                f.write(f"  {sv_type:10}: {file_path}\n")
            f.write("\n")

        # Populations
        populations = demo_data.get('populations', [])
        if populations:
            f.write("\U0001F539 Populations:\n")
            f.write("-" * 40 + "\n")
            for pop in populations:
                f.write(f"  {pop.get('id', 'N/A'):10}: size = {pop.get('initial_size', 'N/A')}\n")
            f.write("\n")

        f.write("\U0001F539 General Configuration:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Output Directory    : {args.output_dir}\n")
        f.write(f"Memory Multiplier   : {config.get('memory_multiplier', 'N/A')}\n\n")

        f.write("\U0001F539 Installed Packages:\n")
        f.write("-" * 40 + "\n")
        for pkg, version in sorted(installed_packages.items(), key=lambda x: x[0].lower()):
            f.write(f"{pkg:30}: {version}\n")
        f.write("\n")

    MScompute(f"Created global simulation recap")


if __name__ == "__main__":
    MScompute(f"Setting up...")
    main()