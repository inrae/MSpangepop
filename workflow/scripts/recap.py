#!/usr/bin/env python3

import os
import yaml
import argparse
import socket
import datetime
import platform
from io_handler import MSsuccess

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
    recap_file = os.path.join(args.output_dir, f"{args.current_run}_params_recap.txt")

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

    # Write the recap file
    with open(recap_file, "w") as f:
        # Header
        f.write(f"🔹 MSpangepop Recap File for: {args.current_run}\n\n")

        # Global Run Information
        f.write("🔹 Global Run Information:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Machine Name       : {machine_name}\n")
        f.write(f"Date/Time          : {current_date}\n")
        f.write(f"Python Version     : {python_version}\n")
        f.write(f"Operating System   : {os_info}\n")
        f.write(f"Current Directory  : {current_dir}\n\n")

        # General Configuration
        f.write("🔹 General Configuration:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Container Registry  : {config.get('container_registry', 'N/A')}\n")
        f.write(f"Output Directory    : {args.output_dir}\n")
        f.write(f"Memory Multiplier   : {config.get('memory_multiplier', 'N/A')}\n")
        f.write(f"SV Type File        : {config.get('sv_type_file', 'N/A')}\n\n")

        # Current Run Parameters
        f.write("🔹 Run Parameters:\n")
        f.write("-" * 30 + "\n")
        samples = config.get("samples", {})
        if args.current_run in samples:
            for param, value in samples[args.current_run].items():
                f.write(f"{param}: {value}\n")
        else:
            f.write(f"⚠️ Warning: No parameters found for {args.current_run}\n")

    MSsuccess(f"Recap file saved: {recap_file}")

if __name__ == "__main__":
    main()
