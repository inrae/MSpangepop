# Warnings / Issues

> **`/!\`:** Act with care; this workflow uses significant memory if you increase the values in `.masterconfig`. We recommend keeping the default settings and running a test first.

> **`/!\`:** For now dont run multiple split at once

> **`/!\`:** The Transduplication and Reciprocal Translocation sections in the `visor_sv_type.yaml` config file are placeholders; do not use them yet.

# How to Use
## A. Running on the CBIB
### 1. Set up
Clone the Git repository and switch to my branch:
```bash
git clone https://forgemia.inra.fr/pangepop/MSpangepop.git
cd MSpangepop
git checkout dev_lpiat
```

### 2. Add your files
- Add a `.fasta.gz` file; an example can be found in the repository.

### 3. Configure the pipeline
- Edit the `.masterconfig` file in the `.config/` directory with your sample information. 
- Edit the `visor_sv_type.yaml` file with the mutations you want.
- Edit line 17 of `job.sh` and line 13 of `./config/snakemake_profile/clusterconfig.yaml` with your email.

### 4. Run the WF
The workflow has two parts: `split` and `simulate`. Always run the split first and once its done (realy quick) run the simulate.
```
sbatch job.sh [split or simulate] dry
```
If no warnings are displayed, run:
```
sbatch job.sh [split or simulate] 
```
> **Nb 1:** to create a visual representation of the workflow, use `dag` instead of `dry`. Open the generated `.dot` file with a [viewer](https://dreampuf.github.io/GraphvizOnline/) that supports the format.

> **Nb 2:** Frist execution of the workflow will be slow since images need to be pulled.

> **Nb 3:** The workflow is in two parts because we want to execute the simulations chromosome by chromosome. Snakemake cannot retrieve the number of chromosomes in one go and needs to index and split first.

## B. Run localy
- Ensure `snakemake` and `singularity` are installed on your machine, then run the workflow:
```
./local_run [split or simulate] dry
```
If the workflow cannot download images from the container registry, install `Docker`, log in with your credentials, and rerun the workflow:
```
docker login -u "<your_username>" -p "<your_token>" "registry.forgemia.inra.fr" 
```

# Workflow
![Dag of the workflow](workflow/dag.svg)

# More informations

The variants generation is inspired by [VISOR](https://github.com/davidebolo1993/VISOR).

You can extract a VCF from the graph using the `vg deconstruct` command. It is not implemented in the pipeline.

# Dependencies
TODO
pandas, msprime, argprase, os, multiprocessing, yaml, Bio.Seq
singularity, snakemake
vg:1.60.0, bcftools:1.12, bgzip:latest, tabix:1.7. 