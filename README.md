# MSpangepop Documentation

MSpangepop is a workflow for simulating variation graphs from coalescent simulations. 

## Documentation Structure

- **[Input Files](doc/input_file.md)** - Requirements and specifications for input data
- **[Master Configuration](doc/configuration.md)** - How to set up configuration files and parameters
- **[Demographic model configuration](doc/create_your_model.md)** - Adapt or create a model
- **[Output Files](doc/output_files.md)** - Description of generated results and outputs
- **[Visualizations](doc/visualisations.md)** - Understanding the generated plots and charts
- **[Advanced Topics](doc/go_even_deeper.md)** - In-depth information for power users

## 🚀 How to Use
### 1. Set up

Clone the Git repository
```bash
git https://forge.inrae.fr/pangepop/MSpangepop 
```

- Create an environement for snakemake (from the provided envfile): 
```bash
conda env create -n wf_env -f .config/wf_env.yaml
```  
> Use Miniforge with the conda-forge channel, see why [here](https://science-ouverte.inrae.fr/fr/offre-service/fiches-pratiques-et-recommandations/quelles-alternatives-aux-fonctionnalites-payantes-danaconda) (french)


### 2. Configure the pipeline for your data

Two elements are needed to run the simulation : 
- The `masterconfig` -> **[Master Configuration](doc/configuration.md)**
- The `demographic_file`  -> **[Demographic model configuration](doc/configuration.md)**

#### To do a quick test : 

Edit the `masterconfig` file in the `.config/` directory with your sample information. -> **[Master Configuration](doc/configuration.md)**

```bash
nano .config/masterconfig.yaml
```

Example config with minimal parameters:
```yaml
samples:
  my_first_run:
    fasta_gz: "small_test_genome.fa.gz"
    chr_n: 1
    demographic_file: "simulation_data/Panmictic_Model.json"
    sv_distribution: {SNP: 50, DEL: 20, INS: 20, INV: 10, DUP: 0}
```
- `fasta_gz` is the input fasta file
- `chr_n` is the number of chromosomes in that file
- `demographic_file` is the demographic scenario the simulation will run on. You can create your own or tailor the ones in `./simulation_data` -> **[Demographic model configuration](doc/configuration.md)**
- `sv_distribution` percentage of each variant type (must sum to 100)



### 3. Run the workflow 
#### On the cluster
- Run the workflow :
```bash
sbatch mspangepop dry # Check for warnings
sbatch mspangepop run # Then
```
> **Nb :** If your account name can't be automatically determined, add it in the `.config/snakemake/profiles/slurm/config.yaml` file.

> **Nb :** Use the command `squeue --format="%.10i %.9P %.6j %.10k %.8u %.2t %.10M %.6D %.20R" -A $user` to see job **names**

#### Localy
```bash
./mspangepop dry # Check for warnings
./mspangepop local-run # Then
```

## ⚙️ Other runing options
```
mspangepop [dry|run|local-run|dag|rulegraph|unlock|touch] [additional snakemake args]
    dry - run in dry-run mode
    run - run the workflow with SLURM
    local-run - run the workflow localy (on a single node)
    dag - generate the directed acyclic graph for the workflow
    rulegraph - generate the rulegraph for the workflow
    unlock - Unlock the directory if snakemake crashed
    touch - Tell snakemake that all files are up to date (use with caution)
    [additional snakemake args] - for any snakemake arg, like --until hifiasm
```

## Support and Development

MSpangepop is developed at INRAe as part of the PangenOak project.


