# Warnings / Issues

> **`/!\`:** Act with care; this workflow uses significant memory if you increase the values in `masterconfig`. We recommend keeping the default settings and running a test first.

# How to Use
## A. Running on a cluster
### 1. Set up
Clone the Git repository
```bash
git clone https://forgemia.inra.fr/pangepop/MSpangepop.git && cd MSpangepop
```

### 2. Add your files
- Add a `.fasta.gz` file; an example can be found in the repository.

### 3. Configure the pipeline
- Edit the `masterconfig` file in the `.config/` directory with your sample information. 
- Edit the `visor_sv_type.yaml` file with the mutations you want.
- Edit `job.sh` with path to the needed modules (`Singularity/Apptainer`, `Miniconda3`)
- Provide the needed conda environement in `job.sh`, under `source activate wf_env`you can create it using :
```bash
conda create -n wf_env -c conda-forge -c bioconda snakemake=8.4.7 snakemake-executor-plugin-slurm
conda init bash
```

### 4. Run the WF
```bash
sbatch job.sh dry
```
If no warnings are displayed, run:
```bash
sbatch job.sh run
```
> **Nb 1:** If the your account name cant be automaticly determined, add it in the `.config/snakemake/profiles/slurm/config.yaml` file.

> **Nb 2:** to create a visual representation of the workflow, use `dag` instead of `dry`. Open the generated `.dot` file with a [viewer](https://dreampuf.github.io/GraphvizOnline/) that supports the format.

## B. Run localy
- Ensure `snakemake` and `Singularity/Apptainer` are installed on your machine, then run the workflow:
```bash
./local_run dry
```
If the workflow cannot download images from the container registry, install `Docker`, log in with your credentials, and rerun the workflow:
```bash
docker login -u "<your_username>" -p "<your_token>" "registry.forgemia.inra.fr" 
```

# More informations
The variants generation is inspired by [VISOR](https://github.com/davidebolo1993/VISOR).


# Dependencies
TODO

Containers :
Miniconda 3, Singularity/Apptainer

Python :
pandas, msprime, Bio.Seq

Workflow :
Snakemake, snakemake-executor-plugin-slurm, Samtool 1.21