# Load the configuration file
configfile: ".config/.masterconfig.yaml"

import os

# Set the output directory for results
output_dir = "results/"

# A default target for Snakemake
rule all:
    input:
        # Ensure that both FAI and FASTA files are included in the all rule
        expand(os.path.join(output_dir, "{sample}_results", "01_split_fai"), sample=config["samples"].keys()) 


# Rule to split the FAI file
rule split_fai:
    input:
        fai=lambda wildcards: config["samples"][wildcards.sample]["fai"]  # Get FAI from config for each sample
    output:
        directory(os.path.join(output_dir, "{sample}_results", "01_split_fai"))  # Directory where split files will be stored
    params:
        out=os.path.join(output_dir, "{sample}_results", "01_split_fai")
    shell:
        """
        mkdir -p {params.out}
        awk '{{print > "{params.out}/" $1 ".fai"}}' {input.fai}
        """

