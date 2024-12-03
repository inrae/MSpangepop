# Load the configuration file
configfile: ".config/masterconfig.yaml"

import os
import yaml

# Retrieve variables from the config file
container_registry = config.get("container_registry", "docker://registry.forgemia.inra.fr/pangepop/mspangepop")
output_dir = config.get("output_dir", "results/")

# Retrieve memory multiplier from config
memory_multiplier = config.get("memory_multiplier", 1)

# Function to load chromosomes dynamically from chr_config.yaml for each sample
def load_chromosomes(sample):
    chr_config_path = os.path.join(output_dir, f"{sample}_results", "01_chromosome_index", "chr_config.yaml")
    if os.path.exists(chr_config_path):
        with open(chr_config_path, 'r') as f:
            chromosomes = yaml.safe_load(f)["chromosomes"]
            return chromosomes
    return []


# Define a function to get the path of the FAI file for each sample and chromosome
def get_fai(wildcards):
    sample = wildcards.sample
    chromosome = wildcards.chromosome
    return os.path.join(output_dir, f"{sample}_results", "02_splited_index", f"{chromosome}.fai")

#### TEST
# Define all possible samples and chromosomes dynamically
samples = list(config["samples"].keys())
chromosomes = []

# Generate chromosomes list for each sample
for sample in samples:
    chromosomes.extend(load_chromosomes(sample))

# Rule to simulate all msprime simulations for all samples and chromosomes
rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.json"),
               sample=samples, chromosome=chromosomes)

rule msprime_simulation:
    input:
        fai=get_fai
    output:
        os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.json")
    params:
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],
        reco_rate=lambda wildcards: config["samples"][wildcards.sample]["recombination_rate"],
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],
        out=lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}_results", "03_msprime_simulation")
    resources:
        mem_mb=lambda wildcards: int(8000 * memory_multiplier),
        time="10:00:00"
    container:
        "/root/MSpangepop/temporary_dependencies/msprime_box.sif"  # Update with the path to the local .sif file
    shell:
        """
        mkdir -p {params.out} &&
        python3 workflow/scripts/tree_generation.py -fai {input.fai} -p {params.pop_size} -m {params.mut_rate} -r {params.reco_rate} -n {params.n} -o {params.out} -c {wildcards.chromosome}
        """
