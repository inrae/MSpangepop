# Load configurations
configfile: ".masterconfig.yaml"  # Load your master configuration file

# Rule to generate all VCFs
rule all:
    input:
        expand("results/{sample}_vcfs/msprime_simulation.vcf", sample=config["samples"].keys()),  # Use keys to expand
        expand("results/{sample}_vcfs/header.vcf", sample=config["samples"].keys())

# Rule to generate VCFs from the Python script
rule generate_tree:
    input:
        fai=lambda wildcards: config["samples"][wildcards.sample]["fai"],  # Get FAI from config for each sample
    output:
        "results/{sample}_vcfs/msprime_simulation.vcf",
        "results/{sample}_vcfs/header.vcf",
    params:
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],  # Get population size from config
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],   # Get mutation rate from config
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],            # Get sample size from config
        outdir="results/{sample}_vcfs",
    shell:
        """
        python3 tree_generation.py -fai {input.fai} -p {params.pop_size} -r {params.mut_rate} -n {params.n} -o {params.outdir}
        """

# Define the samples to run based on the configuration
config["samples"] = {k: v for k, v in config["samples"].items()}  # Ensure samples are directly accessible
