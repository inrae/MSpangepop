# Load the configuration file
configfile: ".config/.masterconfig.yaml"

import os

# Set the output directory for results
output_dir = "results/"

# Rule to define all final outputs
rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "04_final_vcf_generation", "simulated_variants.vcf"),
               sample=config["samples"].keys())

def get_fai(wildcards):
    """
    Retrieve the path of the FAI file for a given sample and chromosome.

    Args:
        wildcards: Wildcards object containing sample and chromosome information.

    Returns:
        str: Path to the corresponding FAI file.
    """
    sample = wildcards.sample
    chromosome = wildcards.chromosome
    return os.path.join(output_dir, f"{sample}_results", "01_split_fai", f"{chromosome}.fai")

rule generate_tree:
    input:
        # Access the correct FAI file using the get_fai function
        fai=get_fai
    output:
        # The expected output VCF file
        os.path.join(output_dir, "{sample}_results", "02_msprime_tree_sim", "{chromosome}_msprime_simulation.vcf")
    params:
        # Use parameters from the config file for population size, mutation rate, and sample size
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],
        out=lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}_results", "02_msprime_tree_sim")
    shell:
        """
        mkdir -p {params.out} &&
        python3 workflow/scripts/tree_generation.py -fai {input.fai} -p {params.pop_size} -r {params.mut_rate} -n {params.n} -o {params.out} -c {wildcards.chromosome}
        """

# Rule to merge all VCF files for each sample
rule merge_vcf:
    input:
        # Collect all chromosome VCFs for a given sample
        vcf_files=expand(os.path.join(output_dir, "{sample}_results", "02_msprime_tree_sim", "{chromosome}_msprime_simulation.vcf"),
                         sample=config["samples"].keys(),
                         chromosome=config["chromosomes"])
    output:
        merged_vcf=os.path.join(output_dir, "{sample}_results", "03_combine_msprime_tree_sim", "combined_msprime_simulation.vcf")
    shell:
        """
        mkdir -p $(dirname {output.merged_vcf}) &&
        bash workflow/scripts/merge_vcf.sh {output.merged_vcf} {input.vcf_files}
        """


rule unzip_fasta:
    input:
        fasta_gz=lambda wildcards: config["samples"][wildcards.sample]["fasta.gz"]
    output:
        temp(output_dir + "{sample}_results/temp/{sample}.fasta")  # Store unzipped FASTA files in a temporary directory
    shell:
        """
        gunzip -c {input.fasta_gz} > {output}
        """


# Rule to generate structural variants using variants_generation.py
rule generate_structural_variants:
    input:
        fai=rules.unzip_fasta.output,
        fasta=output_dir + "{sample}_results/temp/{sample}.fasta",
        vcf=rules.merge_vcf.output,
        yaml=".config/visor_sv_type.yaml"
    output:
        output_dir + "{sample}_results/04_final_vcf_generation/simulated_variants.vcf"
    params:
        outfile=output_dir + "{sample}_results/04_final_vcf_generation/simulated_variants.vcf"
    shell:
        """
        python3 workflow/scripts/variants_generation.py -fai {input.fai} -fa {input.fasta} -v {input.vcf} -y {input.yaml} -o {params.outfile} && 
        rm random_var.tsv
        """
