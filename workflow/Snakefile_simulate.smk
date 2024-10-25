# Load the configuration file
configfile: ".config/.masterconfig.yaml"

import os
import yaml

# Set the output directory for results
output_dir = "results/"

# Function to load chromosomes dynamically from chr_config.yaml for each sample
def load_chromosomes(sample):
    chr_config_path = os.path.join(output_dir, f"{sample}_results", "01_chromosome_index", "chr_config.yaml")
    if os.path.exists(chr_config_path):
        with open(chr_config_path, 'r') as f:
            chromosomes = yaml.safe_load(f)["chromosomes"]
            return chromosomes
    return []

# Rule to define all final outputs
rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "05_final_vcf", "simulated_variants.vcf"),
               sample=config["samples"].keys())

# Define a function to get the path of the FAI file for each sample and chromosome
def get_fai(wildcards):
    sample = wildcards.sample
    chromosome = wildcards.chromosome
    return os.path.join(output_dir, f"{sample}_results", "02_splited_index", f"{chromosome}.fai")

# Rule to generate trees for each chromosome of each sample
rule run_msprime:
    input:
        fai=get_fai
    output:
        os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf")
    params:
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],
        out=lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}_results", "03_msprime_simulation")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/mspangepop_dep:0.0.1"
    shell:
        """
        mkdir -p {params.out} &&
        python3 workflow/scripts/tree_generation.py -fai {input.fai} -p {params.pop_size} -r {params.mut_rate} -n {params.n} -o {params.out} -c {wildcards.chromosome}
        """

# Rule to merge VCF files for each sample by combining VCFs from all chromosomes
rule merge_simulations:
    input:
        vcf_files=lambda wildcards: expand(
            os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf"),
            sample=[wildcards.sample],
            chromosome=load_chromosomes(wildcards.sample)
        )
    output:
        merged_vcf=os.path.join(output_dir, "{sample}_results", "04_merged_msprime_simulation", "combined_msprime_simulation.vcf")
    shell:
        """
        mkdir -p $(dirname {output.merged_vcf}) &&
        bash workflow/scripts/merge_vcf.sh {output.merged_vcf} {input.vcf_files}
        """

# Rule to unzip the FASTA file for each sample
rule unzip_fasta:
    input:
        fasta_gz=lambda wildcards: config["samples"][wildcards.sample]["fasta_gz"]
    output:
        temp(output_dir + "{sample}_results/temp/{sample}.fasta")
    shell:
        """
        gunzip -c {input.fasta_gz} > {output}
        """

# Rule to generate structural variants using variants_generation.py
rule generate_sv:
    input:
        fai=rules.unzip_fasta.output,
        fasta=output_dir + "{sample}_results/temp/{sample}.fasta",
        vcf=rules.merge_simulations.output,
        yaml=".config/visor_sv_type.yaml"
    output:
        output_dir + "{sample}_results/05_final_vcf/simulated_variants.vcf"
    params:
        outfile=output_dir + "{sample}_results/05_final_vcf/simulated_variants.vcf"
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/mspangepop_dep:0.0.1"
    shell:
        """
        python3 workflow/scripts/variants_generation.py -fai {input.fai} -fa {input.fasta} -v {input.vcf} -y {input.yaml} -o {params.outfile} && 
        rm random_var.tsv
        """