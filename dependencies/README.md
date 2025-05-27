This folder contains the dependencies for MSpangepop not yet uploaded to the registery.

# To have the right version of Snakemake

conda create -n wf_env -c conda-forge -c bioconda snakemake=8.4.7

# Run this to update the dependencies box, from the root of the project

rm ./temporary_dependencies/msprime_box.sif && docker build -t msprime_box ./temporary_dependencies/ && singularity build ./temporary_dependencies/msprime_box.sif docker-daemon://msprime_box:latest

# Run test 
singularity exec temporary_dependencies/msprime_box.sif pytest workflow/scripts/test.py 