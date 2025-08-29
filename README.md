# MSpangepop Documentation

MSpangepop is a workflow for simulating variation graphs from coalescent simulations. 

## Documentation Structure

- **[Configuration](doc/configuration.md)** - How to set up configuration files and parameters
- **[Input Files](doc/input_file.md)** - Requirements and specifications for input data
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
- Edit the `masterconfig` file in the `.config/` directory with your sample information. 
```bash
nano .config/masterconfig.yaml
```
- Here you can add the path to your reference genome

Example config : 
```yaml
samples:       
  my_frist_run:                
    fasta_gz: "small_test_genome.fa.gz"  # You can try with this small genome for your first run
    chr_n: 1
    population_size: 5000
    mutation_rate: 1e-5
    recombination_rate: 1e-7
    sample_size: 10
```
#### ⚠️ Important warning :

You can tailor the config file to your dataset with many parameters -> **[Configuration](doc/configuration.md)**

Keep in mind that the parameters of the simulation need to be adjusted to the genome size. 
For large genomes, please start with `1e-11` for `mutation_rate` and `recombination_rate` then go down to your desired value. 

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

## ⚙️ Other running options
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


