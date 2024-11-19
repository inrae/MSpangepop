#!/bin/bash
################################ Slurm options #################################
### prepare_calling_jobs
#SBATCH -J smk_main
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -o slurm_logs/snakemake.%N.%j.out
#SBATCH -e slurm_logs/snakemake.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your.email@here.fr>
################################################################################

# Load the Singularity container runtime
module purge
module load "containers/singularity/3.9.9" || { echo "Failed to load Singularity module"; exit 1; }

# Define paths to your containers
PYTHON_CONTAINER="docker://registry.forgemia.inra.fr/pangepop/mspangepop/python:3.9.7"
SNAKEMAKE_CONTAINER="docker://registry.forgemia.inra.fr/pangepop/mspangepop/snakemake:6.5.1"

# Bind paths
SNG_BIND=$(pwd)
CLUSTER_CONFIG=".config/snakemake_profile/slurm/cluster_config.yml"
MAX_CORES=10
PROFILE=".config/snakemake_profile/slurm"

# Print job details
echo '########################################'
echo 'Job Details:'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job ID:' $SLURM_JOB_ID
echo 'Nodes Allocated:' $SLURM_JOB_NODELIST
echo 'CPUs Per Task:' $SLURM_CPUS_PER_TASK
echo 'Working Directory:' $(pwd)
echo '########################################'


# Run Snakemake inside the Snakemake container
run_snakemake() {
    local snakefile="$1"  # The Snakefile to run
    local option="$2"     # The option for dry run or DAG

    echo "Starting Snakemake workflow with container..."

    singularity exec --bind $SNG_BIND \
        $SNAKEMAKE_CONTAINER \
        snakemake -s "$snakefile" --profile $PROFILE -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" --cluster-config $CLUSTER_CONFIG $option

    # Check if the Snakemake command was successful
    if [ $? -eq 0 ]; then
        echo "$snakefile completed successfully."
    else
        echo "Error: $snakefile failed."
        exit 1
    fi
}

if [ $# -eq 0 ]; then
    echo "Usage: $0 [split|simulate] [dry|dag|run]"
    echo "    split - run the split Snakefile"
    echo "    simulate - run the simulate Snakefile"
    echo "    dry - run the specified Snakefile in dry-run mode"
    echo "    dag - generate DAG for the specified Snakefile"
    echo "    run - run the specified Snakefile normally (default)"
    exit 1
fi

# Determine the workflow and option based on the arguments
workflow="$1"
option="$2"

# Parse options for Snakemake
case "$option" in
    dry)
        snakemake_option="-n -r"
        ;;
    dag)
        snakemake_option="--dag > dag.dot"
        echo "DAG will be generated as dag.png"
        ;;
    *)
        snakemake_option=""
        ;;
esac

case "$workflow" in
    split)
        snakefile="workflow/Snakefile_split.smk"
        ;;
    simulate)
        snakefile="workflow/Snakefile_simulate.smk"
        ;;
    *)
        echo "Invalid workflow: $workflow"
        echo "Usage: $0 [split|simulate] [dry|dag|run]"
        exit 1
        ;;
esac

# Run the Snakemake command
run_snakemake "$snakefile" "$snakemake_option"

# Monitor job queue
squeue -u $USER
