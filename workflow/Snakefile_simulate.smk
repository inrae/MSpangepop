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
        expand(os.path.join(output_dir, "{sample}_results", "06_graph", "graph.gfa"),
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

rule compress_vcf:
    input:
        vcf=os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf")
    output:
        vcf_gz=os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf.gz")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/bgzip:latest"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz}
        """

rule index_vcf:
    input:
        vcf_gz=os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf.gz")
    output:
        index=os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf.gz.tbi")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/tabix:1.7"
    shell:
        """
        tabix -p vcf {input.vcf_gz}
        """

# Rule to merge VCF files for each sample by combining VCFs from all chromosomes
rule merge_simulations:
    input:
        vcf_files=lambda wildcards: expand(
            os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf.gz"),
            sample=[wildcards.sample],
            chromosome=load_chromosomes(wildcards.sample)
        ),
        index_files=lambda wildcards: expand(
            os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "{chromosome}_msprime_simulation.vcf.gz.tbi"),
            sample=[wildcards.sample],
            chromosome=load_chromosomes(wildcards.sample)
        )
    output:
        merged_vcf=os.path.join(output_dir, "{sample}_results", "04_merged_msprime_simulation", "combined_msprime_simulation.vcf.gz")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/bcftools:1.12"
    shell:
        """
        bcftools concat -a -O z -o{output.merged_vcf} -O z {input.vcf_files}
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

rule unzip_vcf:
    input:
        vcf_gz=rules.merge_simulations.output.merged_vcf
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime.vcf")
    shell:
        """
        gunzip -c {input.vcf_gz} > {output}
        """

rule remove_vcf_header:
    input:
        vcf=rules.unzip_vcf.output
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime_no_header.vcf")
    shell:
        """
        awk '!/^#/' {input.vcf} > {output}
        """

# Rule to extract the header from the uncompressed VCF file
rule extract_vcf_header:
    input:
        vcf=rules.unzip_vcf.output
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime_header.vcf")
    shell:
        """
        awk '/^#/' {input.vcf} > {output}
        """

# Rule to generate structural variants using variants_generation.py
rule generate_sv:
    input:
        fai=os.path.join(output_dir, "{sample}_results", "01_chromosome_index", "{sample}_full.fai"),
        fasta=rules.unzip_fasta.output,
        vcf=rules.remove_vcf_header.output,
        yaml=".config/visor_sv_type.yaml"
    output:
        outfile=os.path.join(output_dir, "{sample}_results", "05_final_vcf", "simulated_variants.vcf")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/mspangepop_dep:0.0.1"
    shell:
        """
        python3 workflow/scripts/variants_generation.py -fai {input.fai} -fa {input.fasta} -v {input.vcf} -y {input.yaml} -o {output.outfile} && 
        rm random_var.tsv
        """

# Rule to add the header to the generated VCF output
rule add_header:
    input:
        header=rules.extract_vcf_header.output,
        vcf=rules.generate_sv.output
    output:
        final_output=os.path.join(output_dir, "{sample}_results", "05_final_vcf", "final_simulated_variants.vcf")
    shell:
        """
        cat {input.header} {input.vcf} > {output.final_output}
        """

rule sort_vcf:
    input:
        vcf = rules.add_header.output
    output:
        outfile = os.path.join(output_dir, "{sample}_results", "05_final_vcf", "sorted_simulated_variants.vcf")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/bcftools:1.12"
    shell:
        "bcftools sort {input.vcf} -Oz -o {output}"

rule construct_graph:
    input:
        fasta=rules.unzip_fasta.output, 
        vcf=rules.sort_vcf.output
    output:
        outfile=os.path.join(output_dir, "{sample}_results", "06_graph", "graph.vg")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/vg:1.60.0"
    shell:
        "vg construct -m 2000000000 -r {input.fasta} -v {input.vcf} -f -p > {output}"

rule convert_to_gfa:
    input:
        vg=rules.construct_graph.output
    output:
        outfile=os.path.join(output_dir, "{sample}_results", "06_graph", "graph.gfa")
    container:
        "docker://registry.forgemia.inra.fr/pangepop/mspangepop/vg:1.60.0"
    shell:
        "vg convert -f {input.vg} > {output.outfile}"
