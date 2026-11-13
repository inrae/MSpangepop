# MSpangenome 

MSpangenome is a workflow for simulating pangenome variation graphs from coalescent simulations.\
A simplified description of the algorithm can be found [here](https://lpiat-pages-cfefee.pages-forge.inrae.fr/poster_bims_lucien.pdf).

The official MSpangenome repository can be found at the [INRAE forge](https://forge.inrae.fr/pangepop/MSpangepop).\
A GitHub mirror can be found at [INRAE GitHub](https://github.com/inrae/MSpangepop).\
The mirror is especially useful for people with no Renater account to submit issues.

<img src="doc/pipe.png">

[[_TOC_]]

## Workflow overview

| Stage | Process | Scripts (dev reference) | smk rules (dev reference) |
|-------|---------|---------|--------|
| **1. Setup** | → Validate FASTA/YAML<br>→ Expand configs<br>→ Create index | `input_index.py`<br>`sample_ranges.py`<br>`recap.py` | setup |
| **2. Msprime Simulation** | → Build demographic model<br>→ Run msprime<br>→ Generate visualizations | `msprime_simulation.py`<br>`visualizer_arg.py`<br>`visualizer_tree.py` | msprime_simulation <br>visualization  |
| **3. Preprocessing** | → Split by locus<br>→ Preorder traverse trees<br>→ Assign SV type, length and position | `coalescent_traversal.py`<br>`draw_variants.py`<br>`split_recombination.py` | coalescent_traversal<br>draw_variants<br>split_recombination |
| **4. Graph Creation** | **Initialize:** Build locus ancestral graphs<br>**Mutate:** Apply variants using the MSpangenome library<br>**Save:** Assign IDs → VG unchop each locus → Group loci into subgraph GFAs (100 loci per file by default) | `graph_creation.py`<br>`graph_classes.py`<br>`graph_utils.py`<br>`matrix.py` | graph_creation |
| **5. Merge** | Concatenate the subgraph GFAs → VG unchop the junctions between groups | `merge_subgraphs.py` | gfa_merge |

## Guides and references

|  |  |
|-------|-------------|
| **[Master Configuration](doc/configuration_master.md)** | How to set up configuration files and parameters |
| **[Demographic model](doc/configuration_model.md)** | Adapt or create a model |
| **[Output Files](doc/output_files.md)** | Description of generated results and outputs |
| **[Mutation log](doc/mutation_log.md)** | Format of the per chromosome registry of applied mutations |
| **[Advanced Topics](doc/go_even_deeper.md)** | In-depth information for power users |

---
## Prerequisites

The workflow requires Conda and Apptainer. The following versions have been tested successfully:

- Conda 25.7.0 (Miniforge is recommended)
- Apptainer 1.5.4

### On an HPC cluster

These prerequisites may already be available as environment modules. Check with your sysadmin or consult your cluster documentation.

### On a Linux PC

#### Install Conda

We recommend using Miniforge, a lightweight Conda distribution based on the community-maintained conda-forge channel.

If available, install Miniforge through your Linux distribution's package manager.

If Miniforge is not available or is outdated in your distribution's repositories, follow the installation instructions provided in the official Miniforge documentation:

https://github.com/conda-forge/miniforge

Once installed, check that Conda is available:
`conda --version`

#### Install Apptainer

Install Apptainer using your Linux distribution's package manager, if available. Installation instructions for different Linux distributions are provided in the official Apptainer 1.5 documentation:

https://apptainer.org/docs/admin/1.5/installation.html

Once installed, check that it is available:
`apptainer --version`

### On a macOS or Windows PC

This workflow requires a Linux environment. On macOS, use a supported Linux virtual machine (e.g. [Lima](https://lima-vm.io/)). On Windows, use [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/).

Install Conda and Apptainer inside the Linux environment, following the instructions above.

## How to Use 
### 1. Set up

Clone the Git repository
```bash
git clone https://forge.inrae.fr/pangepop/MSpangepop 
```

- Create an environement for snakemake (from the provided envfile): 
```bash
conda env create -n wf_env -f dependencies/wf_env.yaml
```  

Activate the environment before running the workflow:
```bash
conda activate wf_env
```

### 2. Configure the pipeline for your data

Three elements are needed to run the simulation : 
- The `masterconfig` -> **[Master Configuration](doc/configuration_master.md)**
- The `demographic_file`  -> **[Demographic model configuration](doc/configuration_model.md)**
- A reference FASTA specified by `fasta_gz` (plain or gzip-compressed). The workflow uses the first `chr_n` contigs in the file; one contig is enough for a small test run.

#### To do a quick test : 

Edit the `masterconfig` file in the `.config/` directory with your sample information. (**[Master Configuration](doc/configuration_master.md)**)

```bash
nano .config/masterconfig.yaml
```

Example config:
```yaml
samples:
  test_run:
    model: "simulation_data/Panmictic_Model.json"
    replicates: 1
```
- `model` is the demographic scenario the simulation will run on. You can create your own or tailor the ones in `./simulation_data` (**[Demographic model configuration](doc/configuration_model.md)**)


- **⚠️ Don't want to create your own model?⚠️** Use the provided `Panmictic_Model.json` - simply edit it to specify your genome (`fasta_gz`), then adjust `mutation_rate` and `recombination_rate` (start with low values)




### 3. Run the workflow 
#### On the cluster

- Activate the [Apptainer](https://apptainer.org/docs/user/latest/) / Singularity module from your cluster

- Run the workflow :
```bash
sbatch mspangenome dry # Check for warnings
sbatch mspangenome run # Then
```
> **Nb :** If your account name can't be automatically determined, add it in the `.config/snakemake/profiles/slurm/config.yaml` file.

> **Nb :** Use the command `squeue --format="%.10i %.9P %.6j %.10k %.8u %.2t %.10M %.6D %.20R" -A $user` to see job **names**

> **Nb :** The working directory is mounted in the containers automatically. A genome
> or an output directory living somewhere else needs `--bind`, otherwise the rule that
> reads it fails with a missing file :
> ```bash
> sbatch mspangenome run --bind /work/refs,/scratch/$USER
> ```

#### Locally

Once the [prerequisites](#prerequisites) are installed, run:

```bash
./mspangenome dry # Check for warnings
./mspangenome local-run # Then
```

## Other running options
```
mspangenome [dry|run|local-run|dag|rulegraph|unlock|touch] [-c|--config PATH] [--bind PATHS] [additional snakemake args]
    dry - run in dry-run mode
    run - run the workflow with SLURM
    local-run - run the workflow locally (on a single node)
    dag - generate the directed acyclic graph for the workflow
    rulegraph - generate the rulegraph for the workflow
    unlock - Unlock the directory if snakemake crashed
    touch - Tell snakemake that all files are up to date (use with caution)
    -c, --config PATH - master config file to use (default: .config/masterconfig.yaml)
    --bind PATHS - extra paths to mount in the containers, comma separated, repeatable
    [additional snakemake args] - for any snakemake arg, like --until coalescent_traversal
```

## Support and Development

MSpangenome is developed at INRAE as part of the PangenOak project.


