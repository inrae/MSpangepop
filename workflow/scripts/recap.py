#!/usr/bin/env python3

"""
Author: Lucien Piat
Institution: INRAe
Project: PangenOak

Usage: Create a recap for MSpangepop simulation
"""

import os
import yaml
import argparse
import socket
import datetime
import platform
from io_handler import MSsuccess, MScompute

# Prevent matplotlib from writing to user config dir
os.environ['MPLCONFIGDIR'] = './.config/matplotlib'

# Get versions for all installed packages
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


def main():
    args = parse_args()
    recap_file = os.path.join(args.output_dir, f"{args.current_run}_global_recap.txt")
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    # Global info
    machine_name   = socket.gethostname()
    current_date   = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    python_version = platform.python_version()
    os_info        = f"{platform.system()} {platform.release()}"
    current_dir    = os.getcwd()

    # Gather all installed packages and versions
    installed_packages = {
        dist.project_name: dist.version
        for dist in pkg_resources.working_set
    }

    # Write recap
    with open(recap_file, "w") as f:
        f.write(f"\U0001F539 MSpangepop Recap File for: {args.current_run}\n\n")

        f.write("\U0001F539 Global Run Information:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Machine Name       : {machine_name}\n")
        f.write(f"Date/Time          : {current_date}\n")
        f.write(f"Python Version     : {python_version}\n")
        f.write(f"Operating System   : {os_info}\n")
        f.write(f"Current Directory  : {current_dir}\n\n")

        f.write("\U0001F539 Installed Packages:\n")
        f.write("-" * 30 + "\n")
        for pkg, version in sorted(installed_packages.items(), key=lambda x: x[0].lower()):
            f.write(f"{pkg:30}: {version}\n")
        f.write("\n")

        f.write("\U0001F539 General Configuration:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Output Directory    : {args.output_dir}\n")
        f.write(f"Memory Multiplier   : {config.get('memory_multiplier', 'N/A')}\n")

    MScompute(f"Created global simulation recap")


if __name__ == "__main__":
    MScompute(f"Setting up...")
    main()
