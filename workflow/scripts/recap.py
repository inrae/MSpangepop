#!/usr/bin/env python3

"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage : Create a recap for MSpangepop simulation
"""

import os
import yaml
import argparse
import socket
import datetime
import platform
from io_handler import MSsuccess

# Import necessary libraries to get their versions
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'
import msprime, tskit, pandas, numpy, matplotlib, IPython, yaml, Bio

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a recap file containing simulation parameters."
    )
    parser.add_argument("config_file", help="Path to the configuration YAML file")
    parser.add_argument("output_dir", help="Path to the output directory")
    parser.add_argument("current_run", help="Name of the current run/sample")
    return parser.parse_args()

def main():
    args = parse_args()
    # Define output recap file path
    recap_file = os.path.join(args.output_dir, f"{args.current_run}_global_recap.txt")

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    # Gather global run information
    machine_name   = socket.gethostname()
    current_date   = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    python_version = platform.python_version()
    os_info        = f"{platform.system()} {platform.release()}"
    current_dir    = os.getcwd()

    # Get package versions
    package_versions = {
        "msprime": msprime.__version__,
        "tskit": tskit.__version__,
        "pandas": pandas.__version__,
        "numpy": numpy.__version__,
        "matplotlib": matplotlib.__version__,
        "IPython": IPython.__version__,
        "pyyaml": yaml.__version__,
        "biopython": Bio.__version__,
    }

    # Write the recap file
    with open(recap_file, "w") as f:
        f.write(f"\U0001F539 MSpangepop Recap File for: {args.current_run}\n\n")
        f.write("\U0001F539 Global Run Information:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Machine Name       : {machine_name}\n")
        f.write(f"Date/Time          : {current_date}\n")
        f.write(f"Python Version     : {python_version}\n")
        f.write(f"Operating System   : {os_info}\n")
        f.write(f"Current Directory  : {current_dir}\n\n")
        
        f.write("\U0001F539 Installed Package Versions:\n")
        f.write("-" * 30 + "\n")
        for pkg, version in package_versions.items():
            f.write(f"{pkg:12}: {version}\n")
        f.write("\n")
        
        f.write("\U0001F539 General Configuration:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Container Registry  : {config.get('container_registry', 'N/A')}\n")
        f.write(f"Output Directory    : {args.output_dir}\n")
        f.write(f"Memory Multiplier   : {config.get('memory_multiplier', 'N/A')}\n")
        f.write(f"SV Type File        : {config.get('sv_type_file', 'N/A')}\n\n")
        
        f.write("\U0001F539 Run Parameters:\n")
        f.write("-" * 30 + "\n")

    MSsuccess(f"Created global simulation recap")

if __name__ == "__main__":
    main()
