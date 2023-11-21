configfile: 'config.yaml'

rule all:
	input:
		config["outdir"] + "/msprime.vcf",
		"results/" + config["name"] + ".vcf",
		"results/graph.gfa"

rule fai:
	input:
		fa = config["fa"]
	output:
		fai = config["fa"] + ".fai"
	container:
		"docker://staphb/samtools:1.18"
	shell:
		"samtools faidx {input.fa} --fai-idx {output.fai}"

### msprime coalescent generation
rule msprime:
	input:
		fai = rules.fai.output
	output:
		msvcf = temp(config["outdir"] + "/msprime_simulation.vcf"),
		header = temp(config["outdir"] + "/header.vcf"),
	params:
		population = config["population_size"],
		mut = config["mutation_rate"],
		size = config["sample_size"],
		outdir = config["outdir"]
	container:
		"../mspv.sif"
	shell:
		"python3 /mspv/tree_generation.py -fai {input.fai} -p {params.population} -r {params.mut} -n {params.size} -o {params.outdir}"

rule msprime_vcf:
	input:
		vcf = rules.msprime.output.msvcf,
		header = rules.msprime.output.header
	output:
		config["outdir"] + "/msprime.vcf"
	container:
		"../mspv.sif"
	shell:
		"cat {input.header} {input.vcf} > {output}"

rule bed_vcf:
	input:
		msvcf = rules.msprime.output.msvcf,
		fa = config["fa"],
		fai = rules.fai.output,
		yaml = config["yaml"]
	output:
		vcf = "results/" + config["name"] + ".vcf",
	params:
		out = config["name"]
	container:
		"../mspv_jammy.sif"
	shell:
		"python3 /mspv/variants_generation.py --vcf {input.msvcf} --fasta {input.fa} --fai {input.fai} -y {input.yaml} -o {params.out}"

use rule msprime_vcf as final_vcf with:
	input:
		vcf = rules.bed_vcf.output,
		header = rules.msprime.output.header
	output:
		"results/test_final.vcf"

rule sort_vcf:
	input:
		rules.final_vcf.output
	output:
		"results/test_final_sorted.vcf"
	shell:
		"bcftools sort {input} -Oz -o {output}"

rule graph:
	input:
		vcf = rules.sort_vcf.output,
		fa = config["fa"]
	output:
		graph = "index.giraffe.gbz"
	container:
		"docker://quay.io/vgteam/vg:v1.52.0"
	shell:
		"vg autoindex -r {input.fa} -v {input.vcf} -w giraffe"

rule graph_gfa:
	input:
		rules.graph.output
	output:
		graph = "results/graph.gfa"
	container:
		"docker://quay.io/vgteam/vg:v1.52.0"
	shell:
		"vg convert -f {input} > {output.graph}"

# rule genomes:
# 	input:
# 		graph = rules.output.graph
# 	output:
# 		dir("genomes")
# 	container:
# 		"vg"
# 	shell:
# 		"vg"

# $IMG vg paths -M -x index.giraffe.gbz
# ## convertir les haplotypes en FASTA (donne un multifasta)
# $IMG vg paths -F -x index.giraffe.gbz > paths_index.fa
# ## convertir le graphe au format GFA
# $IMG vg convert -f index.giraffe.gbz > index.giraffe.gfa