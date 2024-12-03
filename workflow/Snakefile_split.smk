# Load the configuration file
configfile: ".config/masterconfig.yaml"

import os
import yaml

# Retrieve variables from the config file
container_registry = config.get("container_registry", "docker://registry.forgemia.inra.fr/pangepop/mspangepop")
output_dir = config.get("output_dir", "results/")

# Retrieve memory multiplier from config
memory_multiplier = config.get("memory_multiplier", 1)

# Rule all ensures final files are generated
rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "02_splited_index"), sample=config["samples"].keys()),
        expand(os.path.join(output_dir, "{sample}_results", "01_chromosome_index", "chr_config.yaml"), sample=config["samples"].keys())

# Rule to generate FAI index for each FASTA file
rule samtools_index_generation:
    input:
        fasta=lambda wildcards: config["samples"][wildcards.sample]["fasta_gz"]
    output:
        fai=os.path.join(output_dir, "{sample}_results", "01_chromosome_index", "{sample}_full.fai")
    threads: 1
    resources:
        mem_mb=lambda wildcards: int(10000 * memory_multiplier),
        time="10:00:00"
    params: 
        out=os.path.join(output_dir, "{sample}_results", "01_chromosome_index")  
    container:
        f"{container_registry}/samtool:1.21"
    shell:
        """
        samtools faidx {input.fasta} &&
        mv {input.fasta}.fai {output.fai} && echo 'indexing successful' &&
        rm {input.fasta}.gzi || true 
        """

# Rule to split the FAI file into separate chromosome files
rule awk_index_split:
    input:
        fai=rules.samtools_index_generation.output.fai
    output:
        directory(os.path.join(output_dir, "{sample}_results", "02_splited_index"))
    resources:
        mem_mb=lambda wildcards: int(5000 * memory_multiplier),
        time="10:00:00"
    params:
        out=os.path.join(output_dir, "{sample}_results", "02_splited_index")
    shell:
        """
        mkdir -p {params.out}
        awk '{{print > "{params.out}/" $1 ".fai"}}' {input.fai}
        """

# Rule to create chromosome configuration YAML with chromosome names
rule yaml_index_conversion:
    input:
        fai=rules.samtools_index_generation.output.fai
    output:
        yaml=os.path.join(output_dir, "{sample}_results", "01_chromosome_index", "chr_config.yaml")
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="10:00:00"
    shell:
        """
        awk '{{print $1}}' {input.fai} | \
        awk 'BEGIN {{print "chromosomes:"}} {{print "  - \\"" $1 "\\""}}' > {output.yaml}
        """