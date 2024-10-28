# Warnings

> **`/!\`:** Act with care; this workflow uses significant memory if you increase the values in `.masterconfig`. We recommend keeping the default settings and running a test first.

> **`/!\`:** For now workflow only tested with SPN generation

> **`/!\`:** For now dont run multiple split at once


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
- Add a `.fasta.gz` file; examples can be found in the repository

### 3. Configure the pipeline
- Edit the `.masterconfig` file and the `visor_sv_type.yaml` in the `.config/` directory to suit your needs. 
- Edit the `job.sh` and the `./config/snakemake_profile/clusterconfig.yaml` with your email.

### 4. Run the WF
The workflow has two parts: `split` and `simulate`. Always run the split first and once its done (realy quick) run the simulate. To do so run the following commands:
```
sbatch job.sh [split, simulate] [dry, dag]
```
If no warnings are displayed, run:
```
sbatch job.sh [split, simulate] 
```
> **Nb 1:** to create a visual representation of the workflow, use [dag]. Open the generated `.dot` file with a [viewer](https://dreampuf.github.io/GraphvizOnline/) that supports the format.

> **Nb 2:** Frist execution of the workflow will be slow since images need to be pulled.
## B. Run localy
- Ensure `snakemake` and `singularity` are installed on your machine.
- Modify the `.masterconfig` file and `visor_sv_type.yaml` in the `.config/` directory as needed.

```
./local_run [split, simulate] [dry, dag]
```

If the workflow cannot download images from the container registry, install `Docker`, log in with your credentials, and rerun the workflow:
```
docker login -u "<your_username>" -p "<your_token>" "registry.forgemia.inra.fr" 
```

# Workflow
![Dag of the workflow](workflow/dag.svg)


# More informations

The variants generation is inspired by [VISOR](https://github.com/davidebolo1993/VISOR).

Each variant type has a size distribution file (bins = 100 bp) in folder `sv_distributions`. The data was extracted from [An integrated map of structural variation in 2,504 human genomes (Sudmant, et al. 2015)](https://www.nature.com/articles/nature15394). The distributions are used to randomly sample each structural variant size.

You can extract a VCF from the graph using the `vg deconstruct` command. It is not implemented in the pipeline.

# Dependencies
TODO
pandas, msprime, argprase, os, multiprocessing, yaml, Bio.Seq
singularity, snakemake
vg:1.60.0, bcftools:1.12, bgzip:latest, tabix:1.7. 