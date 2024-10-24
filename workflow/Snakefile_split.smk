# Load the configuration file
configfile: ".config/.masterconfig.yaml"

import os

# Set the output directory for results
output_dir = "results/"


rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "02_split_fai"), sample=config["samples"].keys()),
        expand(os.path.join(output_dir, "{sample}_results", "01_all_chromosomes", "chr_config.yaml"), sample=config["samples"].keys()) 

rule generate_fai:
    input:
        fasta=lambda wildcards: config["samples"][wildcards.sample]["fasta_gz"]
    output:
        fai=os.path.join(output_dir, "{sample}_results", "01_all_chromosomes", "{sample}_full.fai")  
    params: 
        out=os.path.join(output_dir, "{sample}_results", "01_all_chromosomes")  
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/samtool:1.21"
    shell:
        """
        samtools faidx {input.fasta} &&
        mv {input.fasta}.fai {output.fai} &&
        rm {input.fasta}.gzi || true
        """

rule split_fai:
    input:
        fai=rules.generate_fai.output.fai
    output:
        directory(os.path.join(output_dir, "{sample}_results", "02_split_fai"))
    params:
        out=os.path.join(output_dir, "{sample}_results", "02_split_fai")
    shell:
        """
        mkdir -p {params.out}
        awk '{{print > "{params.out}/" $1 ".fai"}}' {input.fai}
        """

rule create_chr_config:
    input:
        fai=rules.generate_fai.output.fai
    output:
        yaml=os.path.join(output_dir, "{sample}_results", "01_all_chromosomes", "chr_config.yaml")
    shell:
        """
        bash workflow/scripts/fai2yaml.sh {input.fai} {output.yaml}
        """