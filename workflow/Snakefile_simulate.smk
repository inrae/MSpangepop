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

# Rule to define all final outputs
rule all:
    input:
        expand(os.path.join(output_dir, "{sample}_results", "04_generated_variants", "{sample}_simulated_variants.vcf.gz"),
               sample=config["samples"].keys()) + 
        expand(os.path.join(output_dir, "{sample}_results", "05_vg_graph", "{sample}_vg_graph.gfa"),
               sample=config["samples"].keys()) + 
        expand(os.path.join(output_dir, "{sample}_results", "06_graph_paths", "{sample}_paths.fasta.gz"),
               sample=config["samples"].keys())

# Define a function to get the path of the FAI file for each sample and chromosome
def get_fai(wildcards):
    sample = wildcards.sample
    chromosome = wildcards.chromosome
    return os.path.join(output_dir, f"{sample}_results", "02_splited_index", f"{chromosome}.fai")

# Rule to run msprime simulation
rule msprime_simulation:
    input:
        fai=get_fai
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "{chromosome}_msprime_simulation.vcf"))
    params:
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],
        reco_rate=lambda wildcards: config["samples"][wildcards.sample]["recombination_rate"],
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],
        out=lambda wildcards: os.path.join(output_dir, f"{wildcards.sample}_results", "temp")
    resources:
        mem_mb=lambda wildcards: int(8000 * memory_multiplier),
        time="10:00:00"
    container:
        f"{container_registry}/mspangepop_dep:0.0.1"
    shell:
        """
        mkdir -p {params.out} &&
        python3 workflow/scripts/tree_generation.py -fai {input.fai} -p {params.pop_size} -m {params.mut_rate} -r {params.reco_rate} -n {params.n} -o {params.out} -c {wildcards.chromosome}
        """

# Rule to compress the msprime simulation VCF file
rule bgzip_compress_simvcf:
    input:
        vcf=rules.msprime_simulation.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "{chromosome}_msprime_simulation.vcf.gz"))
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="01:00:00"
    container:
        f"{container_registry}/bgzip:latest"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        """

# Rule to index the compressed VCF file
rule tabix_index_simvcf:
    input:
        vcf_gz=rules.bgzip_compress_simvcf.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "{chromosome}_msprime_simulation.vcf.gz.tbi"))
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:30:00"
    container:
        f"{container_registry}/tabix:1.7"
    shell:
        """
        tabix -p vcf {input.vcf_gz}
        """


# Rule to merge VCF files for each sample by combining VCFs from all chromosomes
rule bcftools_concat_simvcfs:
    input:
        vcf_files=lambda wildcards: expand(
            os.path.join(output_dir, "{sample}_results", "temp", "{chromosome}_msprime_simulation.vcf.gz"),
            sample=[wildcards.sample],
            chromosome=load_chromosomes(wildcards.sample)
        ),
        index_files=lambda wildcards: expand(
            os.path.join(output_dir, "{sample}_results", "temp", "{chromosome}_msprime_simulation.vcf.gz.tbi"),
            sample=[wildcards.sample],
            chromosome=load_chromosomes(wildcards.sample)
        )
    output:
        merged_vcf=os.path.join(output_dir, "{sample}_results", "03_msprime_simulation", "msprime_simulation.vcf.gz")
    resources:
        mem_mb=lambda wildcards: int(5000 * memory_multiplier),
        time="01:30:00"
    container:
        f"{container_registry}/bcftools:1.12"
    shell:
        """
        bcftools concat -a -O z -o {output.merged_vcf} {input.vcf_files}
        """

# Rule to unzip the FASTA file for each sample
rule gunzip_fasta:
    input:
        fasta_gz=lambda wildcards: config["samples"][wildcards.sample]["fasta_gz"]
    output:
        temp(output_dir + "{sample}_results/temp/{sample}.fasta")
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:30:00"
    shell:
        """
        gunzip -c {input.fasta_gz} > {output}
        """

# Rule to unzip the simulated VCF
rule gunzip_simvcf:
    input:
        vcf_gz=rules.bcftools_concat_simvcfs.output.merged_vcf
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime.vcf")
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:30:00"
    shell:
        """
        gunzip -c {input.vcf_gz} > {output}
        """

# Rule to remove the header from the VCF file
rule awk_header_removal:
    input:
        vcf=rules.gunzip_simvcf.output
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime_no_header.vcf")
    resources:
        mem_mb=lambda wildcards: int(500 * memory_multiplier),
        time="00:15:00"
    shell:
        """
        awk '!/^#/' {input.vcf} > {output}
        """

# Rule to extract the header from the uncompressed VCF file
rule awk_header_extraction:
    input:
        vcf=rules.gunzip_simvcf.output
    output:
        temp(output_dir + "{sample}_results/temp/{sample}_msprime_header.vcf")
    resources:
        mem_mb=lambda wildcards: int(500 * memory_multiplier),
        time="00:15:00"
    shell:
        """
        awk '/^#/' {input.vcf} > {output}
        """

# Rule to generate structural variants using the variants_generation.py script
rule generate_variants:
    input:
        fai=os.path.join(output_dir, "{sample}_results", "01_chromosome_index", "{sample}_full.fai"),
        fasta=rules.gunzip_fasta.output,
        vcf=rules.awk_header_removal.output,
        yaml=".config/visor_sv_type.yaml"
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "simulated_variants.vcf"))
    resources:
        mem_mb=lambda wildcards: int(8000 * memory_multiplier),
        time="02:00:00"
    container:
        f"{container_registry}/mspangepop_dep:0.0.1"
    shell:
        """
        python3 workflow/scripts/generate_variant.py --fai {input.fai} --fasta {input.fasta} --vcf {input.vcf} --yaml {input.yaml} --output {output}
        """

# Rule to sort the VCF header
rule sort_header:
    input:
        header=rules.awk_header_extraction.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "sorted_header.vcf"))
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:20:00"
    container:
        f"{container_registry}/mspangepop_dep:0.0.1"
    shell:
        """
        python3 workflow/scripts/sort_vcf_header.py {input.header} {output}
        """

# Rule to add the header to the VCF output
rule cat_add_header:
    input:
        header=rules.sort_header.output,
        vcf=rules.generate_variants.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "final_simulated_variants.vcf"))
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:20:00"
    shell:
        """
        cat {input.header} {input.vcf} > {output}
        """

# Rule to sort and compress the final VCF
rule bcftools_sort_finalvcf:
    input:
        vcf=rules.cat_add_header.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "sorted_simulated_variants.vcf"))
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="00:30:00"
    container:
        f"{container_registry}/bcftools:1.12"
    shell:
        """
        bcftools sort {input.vcf} -Oz -o {output}
        """

# Rule to construct the graph
rule vg_construct_graph:
    input:
        fasta=rules.gunzip_fasta.output,
        vcf=rules.bcftools_sort_finalvcf.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "05_vg_graph", "graph.vg"))
    resources:
        mem_mb=lambda wildcards: int(8000 * memory_multiplier),
        time="02:00:00"
    container:
        f"{container_registry}/vg:1.60.0"
    params:
        memory = lambda wildcards: int(7000 * memory_multiplier),
    shell:
        """
        vg construct -m {params.memory} -r {input.fasta} -v {input.vcf} -f -p > {output}
        """

# Rule to convert VG graph to GFA format
rule vg_convert_to_gfa:
    input:
        vg=rules.vg_construct_graph.output
    output:
        outfile=os.path.join(output_dir, "{sample}_results", "05_vg_graph", "{sample}_vg_graph.gfa")
    resources:
        mem_mb=lambda wildcards: int(4000 * memory_multiplier),
        time="01:00:00"
    container:
        f"{container_registry}/vg:1.60.0"
    shell:
        """
        vg convert -f {input.vg} > {output.outfile}
        """

# Rule to compress VCF for Giraffe
rule bgzip_compress_finalvcf:
    input:
        vcf=rules.bcftools_sort_finalvcf.output
    output:
        os.path.join(output_dir, "{sample}_results", "04_generated_variants", "{sample}_simulated_variants.vcf.gz")
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="00:30:00"
    container:
        f"{container_registry}/bgzip:latest"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        """

# Rule to create Giraffe index
rule vg_autoindex_giraffe:
    input:
        fasta=rules.gunzip_fasta.output,
        vcf_gz=rules.bgzip_compress_finalvcf.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "index.giraffe.gbz"))
    params:
        out=os.path.join(output_dir, "{sample}_results", "temp", "index")
    resources:
        mem_mb=lambda wildcards: int(12000 * memory_multiplier),
        time="03:00:00"
    container:
        f"{container_registry}/vg:1.60.0"
    shell:
        """
        vg autoindex -r {input.fasta} -v {input.vcf_gz} -w giraffe -p {params.out}
        """

# Rule to convert GBZ to GFA format
rule vg_convert_gbz_to_gfa:
    input:
        gbz=rules.vg_autoindex_giraffe.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "giraffe_graph.gfa"))
    resources:
        mem_mb=lambda wildcards: int(4000 * memory_multiplier),
        time="01:30:00"
    container:
        f"{container_registry}/vg:1.60.0"
    shell:
        """
        vg convert -f {input.gbz} > {output}
        """

# Rule to extract paths to a FASTA file
rule vg_paths:
    input:
        gbz=rules.vg_convert_gbz_to_gfa.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "{sample}_paths_dirty.fasta"))
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="01:00:00"
    container:
        f"{container_registry}/vg:1.60.0"
    shell:
        """
        vg paths -F -x {input.gbz} > {output}
        """

# Rule to remove the reference paths
rule remove_reference:
    input:
        fasta=rules.vg_paths.output
    output:
        temp(os.path.join(output_dir, "{sample}_results", "temp", "{sample}_paths_.fasta"))
    resources:
        mem_mb=lambda wildcards: int(1000 * memory_multiplier),
        time="00:20:00"
    shell:
        """
        ./workflow/scripts/remove_reference.sh {input.fasta} {output}
        """

# Rule to compress the filtered FASTA file
rule bgzip_compress_finalfasta:
    input:
        fasta=rules.remove_reference.output
    output:
        os.path.join(output_dir, "{sample}_results", "06_graph_paths", "{sample}_paths.fasta.gz")
    resources:
        mem_mb=lambda wildcards: int(2000 * memory_multiplier),
        time="00:30:00"
    container:
        f"{container_registry}/bgzip:latest"
    shell:
        """
        bgzip -c {input.fasta} > {output}
        """