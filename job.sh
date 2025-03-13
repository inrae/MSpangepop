#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH -o slurm_logs/out_job_%j.out
#SBATCH -e slurm_logs/err_job_%j.err
#SBATCH --time=80:00:00
#SBATCH -J MSpangepop
#SBATCH --mem=10G

# Function to display usage instructions
usage() {
    echo "Usage: $0 [split|simulate] [dry|dag|run]"
    echo "  [split|simulate]: The workflow you want to run"
    echo "  [dry|dag|run]: The option to execute (dry-run, generate DAG, or actual run)"
    echo "Examples:"
    echo "  $0 split dry     # Run the split workflow in dry-run mode"
    echo "  $0 simulate run  # Run the simulate workflow"
    exit 1
}

# Update this with the path to your images
echo 'Loading modules'
module purge
module load containers/Apptainer/1.2.5 
module load devel/Miniconda/Miniconda3

echo 'Activating environment'
source activate wf_env

echo 'Starting Snakemake workflow'

run_snakemake() {
    local snakefile="$1"  
    local option="$2"     
    echo "Starting $snakefile..."

    # Execute the Snakemake command with the specified option
    if [[ "$option" == "dry" ]]; then
        snakemake -s "$snakefile" -c $(nproc) --dry-run
    elif [[ "$option" == "dag" ]]; then
        snakemake -s "$snakefile" -c $(nproc) --dag > workflow.dot
        echo "DAG has been generated as workflow.svg"
        return
    elif [[ "$option" == "run" || -z "$option" ]]; then
        snakemake -s "$snakefile" --workflow-profile ./.config/snakemake/profiles/slurm
    else
        echo "Invalid option: $option"
        usage  # Display usage if option is invalid
    fi


}

# Determine the workflow and option based on the arguments
workflow="$1"
option="$2"

# Validate arguments and call the usage function if invalid
if [[ -z "$workflow" ]]; then
    echo "Error: Missing required workflow argument."
    usage  # Display usage if workflow is missing
fi

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
        usage  # Display usage if workflow is invalid
        ;;
esac

# Run the specified Snakefile with the provided option
run_snakemake "$snakefile" "$option"
