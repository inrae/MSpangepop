#!/bin/bash
################################ Slurm options #################################
### prepare_calling_jobs
#SBATCH -J smk_main
### Max run time "hours:minutes:seconds"
#SBATCH --time=96:00:00
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
#SBATCH --mail-user=lucien.piat@inare.fr
################################################################################

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job ID:' $SLURM_JOB_ID
echo 'Number of nodes assigned to job:' $SLURM_JOB_NUM_NODES
echo 'Total number of cores for job (?):' $SLURM_NTASKS
echo 'Number of requested cores per node:' $SLURM_NTASKS_PER_NODE
echo 'Nodes assigned to job:' $SLURM_JOB_NODELIST
echo 'Number of CPUs assigned for each task:' $SLURM_CPUS_PER_TASK
echo 'Directory:' $(pwd)
# Detail Information:
echo 'scontrol show job:'
scontrol show job $SLURM_JOB_ID
echo '########################################'

# Function to load modules
load_modules() {
    module purge  # Clear any previously loaded modules

    # Loop through each module and load it
    for module_name in "$@"; do
        module load "$module_name"
    done
}

# Here specify the modules to load and their path
load_modules "python/3.9.7" "snakemake/6.5.1" 

### variables
SNG_BIND="/mnt/cbib/pangenoak_trials/MSpangepop/"
CLUSTER_CONFIG=".config/snakemake_profile/slurm/cluster_config.yml"
MAX_CORES=10
PROFILE=".config/snakemake_profile/slurm"

echo 'Starting Snakemake workflow'


run_snakemake() {
    local snakefile="$1"  # The Snakefile to run
    local option="$2"     # The option for dry run or DAG

    echo "Starting $snakefile..."

    # Execute the Snakemake command with the specified option
    if [[ "$option" == "dry" ]]; then
        snakemake -s "$snakefile" --profile $PROFILE -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" --cluster-config $CLUSTER_CONFIG -n -r
    elif [[ "$option" == "dag" ]]; then
        snakemake -s "$snakefile" --profile $PROFILE -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" --cluster-config $CLUSTER_CONFIG --dag > dag.dot
        echo "DAG has been generated as dag.png"
        return
    else
        snakemake -s "$snakefile" --profile $PROFILE -j $MAX_CORES --use-singularity --singularity-args "-B $SNG_BIND" --cluster-config $CLUSTER_CONFIG
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