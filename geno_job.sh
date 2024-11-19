#!/bin/bash
################################ Slurm options #################################
### prepare_calling_jobs
#SBATCH -J smk_main
### Max run time "hours:minutes:seconds"
#SBATCH --time=95:00:00
#SBATCH --ntasks=1 #nb of processes
#SBATCH --cpus-per-task=1 # nb of cores for each process(1 process)
#SBATCH --mem=10G # max of memory (-m) 
### Requirements nodes/servers (default: 1)
#SBATCH --nodes=1
### Requirements cpu/core/task (default: 1)
#SBATCH --ntasks-per-node=1
#SBATCH -o slurm_logs/snakemake.%N.%j.out
#SBATCH -e slurm_logs/snakemake.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<your.email@here.fr>
################################################################################

# This script is used as temporary sollution for runing MSpangepop on the genotoul HPC

# Function to load modules
load_modules() {
    module purge  # Clear any previously loaded modules

    # Loop through each module and load it
    for module_name in "$@"; do
        module load "$module_name"
    done
}

load_modules "containers/singularity/3.9.9" "bioinfo/Snakemake/8.3.1"

SNG_BIND=$(pwd)
MAX_CORES=10

echo 'Starting Snakemake workflow'

run_snakemake() {
    local snakefile="$1"  # The Snakefile to run
    local option="$2"     # The option for dry run or DAG

    echo "Starting $snakefile..."

    # Execute the Snakemake command with the specified option
    if [[ "$option" == "dry" ]]; then
        snakemake -s "$snakefile" -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" -n
    elif [[ "$option" == "dag" ]]; then
        snakemake -s "$snakefile" -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" --dag > dag.dot
        echo "DAG has been generated as dag.png"
        return
    else
        snakemake -s "$snakefile" -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND"
    fi

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
    echo ""
    echo "This script is used as temporary sollution for runing MSpangepop on the genotoul HPC"
    echo "Please use job.sh on any other HPC"
    exit 1
fi

# Determine the workflow and option based on the arguments
workflow="$1"
option="$2"

# Run the specified Snakefile based on user input
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

# Run the specified Snakefile with the provided option
run_snakemake "$snakefile" "$option"