# Load the configuration file
configfile: ".config/.masterconfig.yaml"

import os

# Set the output directory for results
output_dir = "results/"


rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "02_split_fai"), sample=config["samples"].keys()) 

rule generate_fai:
    input:
        fasta=lambda wildcards: config["samples"][wildcards.sample]["fasta_gz"]  # Get FASTA file from config
    output:
        fai=os.path.join(output_dir, "{sample}_results", "01_full_fai", "{sample}_full.fai")  
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/samtool:1.21"
    shell:
        """
        samtools faidx {input.fasta}
        mv {input.fasta}.fai {output.fai}
        """

# Rule to split the FAI file
rule split_fai:
    input:
        fai=rules.generate_fai.output.fai
    output:
        directory(os.path.join(output_dir, "{sample}_results", "02_split_fai"))  # Directory where split files will be stored
    params:
        out=os.path.join(output_dir, "{sample}_results", "02_split_fai")
    shell:
        """
        mkdir -p {params.out}
        awk '{{print > "{params.out}/" $1 ".fai"}}' {input.fai}
        """

