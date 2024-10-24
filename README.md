# How to Use

### 1. Set Up

You can run the workflow directly after downloading; no need for additional libraries; everything is included in a Singularity image.

```bash
git clone https://forgemia.inra.fr/pangepop/MSpangepop.git
cd MSpangepop
git checkout dev_lpiat
```
Clone the branch

### 2. Add your files
- Add a `.fasta.gz` file
- Add a `.fai` file 

you will find examples in the repo. 

### 3. Set up the pipeline
- Edit the `.masterconfig` file and the `visor_sv_type.yaml` in the `.config/` directory to suit your needs. 
- Edit the `job.sh` and the `./config/snakemake_profile/clusterconfig.yaml` with your email.

### 4. Run the WF
```
sbatch job.sh [split, simulate] [dry, dag]
```
If no warnings are displayed
```
sbatch job.sh [split, simulate] 
```

Nb: you can also run localy, use ./local_run [split, simulate] [dry, dag] and make sure you have snakemake and singularity instaled on the local machine

/!\ Act with care; this workflow is a proper memory hog if you increase the values in the .masterconfig too much.

### Help ! It is complaning about the singularity image 
If you want to run localy, make sure you are loged in the forgemia registery before executing the pipeline.
```
docker login -u "<your_username>" -p "<your_token>" "registry.forgemia.inra.fr" 
```

# More informations

The variants generation is inspired by [VISOR](https://github.com/davidebolo1993/VISOR). YAML template is available in `VISOR_random_bed` folder. Modify the YAML input to set the percentage of each variant you want to simulate (must equal 100).

Each variant type has a size distribution file (bins = 100 bp) in folder `sv_distributions`. The data was extracted from [An integrated map of structural variation in 2,504 human genomes (Sudmant, et al. 2015)](https://www.nature.com/articles/nature15394).

The distributions are used to randomly sample each structural variant size.

## Create exact data with vg
In the `vg_extact_data` folder.

Snakemake/Singularity pipeline to get a pangenome in GFA format and a FASTA with all individuals from the VCF. Starts from a reference FASTA and a VCF to specify in the `config.yaml` file, with a name for the output (**WARNING**: think to set the correct email address in `config.yaml`).

Create directory or modify the SBATCH options in `job.sh` (**WARNING**: think to set the correct email address in `job.sh` if you want to receive the slurm emails).
```
mkdir -p slurm_logs
```
On SLURM cluster, run `sbatch job.sh dry` for a dry run or `sbatch job.sh` directly. Adjust the `SNG_BIND` variable if files are not found and the snakemake profile as necessary for performance.

You can extract a VCF from the graph using the `vg deconstruct` command. It is not implemented in the pipeline.
