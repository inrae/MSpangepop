#!/bin/bash

## TMP config to run on the CBIB
#SBATCH --job-name=asm4pg
#SBATCH --ntasks=25
#SBATCH --mem=300G
#SBATCH -o slurm_logs/out_job_%j.out
#SBATCH -e slurm_logs/err_job_%j.err

# Written by Lucien Piat at INRAe
# Use this script to run MSpangepop localy or on a single HPC node
# 03/02/25

SNG_BIND=$(pwd)

run_snakemake() {
    local option="$1"

    case "$option" in
        dry)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) -n 
            ;;
        dag)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) --dag > dag.dot
            if [ $? -eq 0 ]; then
                echo "✅ MSpangepop -> DAG has been successfully generated as dag.dot"
            else
                echo "❌ MSpangepop -> Error: Failed to generate DAG."
                exit 1
            fi
            ;;
        rulegraph)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) --rulegraph > rulegraph.dot
            if [ $? -eq 0 ]; then
                echo "✅ MSpangepop -> Rulegraph has been successfully generated as rulegraph.dot"
            else
                echo "❌ MSpangepop -> Error: Failed to generate Rulegraph."
                exit 1
            fi
            ;;
        unlock)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) --unlock
            ;;
        run)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) -k -q rules
            ;;
        verbose)
            snakemake --use-singularity --singularity-args "-B $SNG_BIND" -j $(nproc) -k --verbose
            ;;
        *)
            echo "Invalid option: $option"
            echo "Usage: $0 [dry|run|dag|rulegraph|unlock|verbose]"
            exit 1
            ;;
    esac
}

# Verify arguments
if [ $# -ne 1 ] || [ "$1" == "help" ]; then
    echo "Use this script to run MSpangepop localy or on a single HPC node"
    echo ""
    echo "Usage: $0 [dry|run|dag|rulegraph|unlock]"
    echo "    dry - run in dry-run mode"
    echo "    run - run normally"
    echo "    dag - generate the directed acyclic graph"
    echo "    rulegraph - generate the rulegraph"
    echo "    unlock - unlock the directory if snakemake crashed"
    echo "    verbose - run in verbose debug mode"
    exit 1
fi

# Execute the function with the provided option
echo "🔹 MSpangepop -> Starting snakemake workflow"
run_snakemake "$1"
