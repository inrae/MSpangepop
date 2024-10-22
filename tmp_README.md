# How to Use

**Important:** Make sure you are not on an HPC (High-Performance Computing) system, as this workflow is not yet designed for SLURM.

### 1. Set Up

Ensure you have the following packages installed: `msprime`, `pandas`, `defopt`, `numpy`, `biopython`, `pyyaml`, `snakemake`, and `gzip`.

```bash
sudo apt update -y
sudo apt install python3 -y
python3 --version
sudo apt install python3-pip -y
pip install msprime pandas defopt numpy biopython pyyaml
sudo apt install gzip
```
Clone the branch

```bash
git clone -b dev_lpiat https://forgemia.inra.fr/pangepop/MSpangepop.git
cd MSpangepop/
```
### 2. Add your files
- Add a `.fasta.gz` file
- Add a `.fai` file 

you will find examples in the repo. 

### 3. Set up the pipeline
Edit the `.masterconfig` file and the `visor_sv_type.yaml` in the `.config/` directory to suit your needs.

### 4. Run the WF
```
snakemake -c1 -n
```
If no warnings are displayed
```
snakemake -c1
```
/!\ Act with care; this workflow is a proper memory hog if you increase the values in the .masterconfig too much.