#!/bin/bash
#Script to run localy the workflow DO NOT USE AS IS ON A CLUSTER!

SNG_BIND=$(pwd)

run_snakemake() {
    local snakefile="$1"  # The Snakefile to run
    local option="$2"     # The option for dry run or DAG

    echo "Starting $snakefile..."

    # Execute the Snakemake command with the specified option
    if [[ "$option" == "dry" ]]; then
        snakemake -s "$snakefile" --use-singularity --singularity-args "-B $SNG_BIND" -c4 -n
    elif [[ "$option" == "dag" ]]; then
                snakemake -s "$snakefile" --use-singularity --singularity-args "-B $SNG_BIND" -c4 --dag > dag.dot
        echo "DAG has been generated as dag.png"
        return
    else
        snakemake -s "$snakefile" --use-singularity --singularity-args "-B $SNG_BIND" -c4
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


