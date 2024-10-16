configfile: ".masterconfig.yaml"  

config["samples"] = {k: v for k, v in config["samples"].items()}

output_dir = "results/"

# Final rule to generate all structural variants
rule all:
    input:
        expand(output_dir + "{sample}_results/02_variant_generation/{sample}_variants.vcf", sample=config["samples"].keys())

# Rule to unzip fasta.gz files to a temporary folder
rule unzip_fasta:
    input:
        fasta_gz=lambda wildcards: config["samples"][wildcards.sample]["fasta.gz"]
    output:
        temp(output_dir + "{sample}_results/temp/{sample}.fa")  # Store unzipped FASTA files in a temporary directory
    shell:
        """
        gunzip -c {input.fasta_gz} > {output}
        """

# Rule to generate VCFs from the Python script
rule generate_tree:
    input:
        fai=lambda wildcards: config["samples"][wildcards.sample]["fai"],  # Get FAI from config for each sample
    output:
        output_dir + "{sample}_results/01_genet_tree/msprime_simulation.vcf"
    params:
        pop_size=lambda wildcards: config["samples"][wildcards.sample]["population_size"],  # Get population size from config
        mut_rate=lambda wildcards: config["samples"][wildcards.sample]["mutation_rate"],   # Get mutation rate from config
        n=lambda wildcards: config["samples"][wildcards.sample]["sample_size"],            # Get sample size from config
        out = output_dir + "{sample}_results/01_genet_tree"
    shell:
        """
        python3 tree_generation.py -fai {input.fai} -p {params.pop_size} -r {params.mut_rate} -n {params.n} -o {params.out}
        """

# Rule to generate structural variants using variants_generation.py
rule generate_structural_variants:
    input:
        fai=lambda wildcards: config["samples"][wildcards.sample]["fai"],
        fasta=output_dir + "{sample}_results/temp/{sample}.fa",
        vcf=output_dir + "{sample}_results/01_genet_tree/msprime_simulation.vcf",
        yaml="visor_sv_type.yaml"
    output:
        output_dir + "{sample}_results/02_variant_generation/{sample}_variants.vcf"
    params:
        outfile=output_dir + "{sample}_results/02_variant_generation/{sample}_variants.vcf"
    shell:
        """
        python3 variants_generation.py -fai {input.fai} -fa {input.fasta} -v {input.vcf} -y {input.yaml} -o {params.outfile}
        """