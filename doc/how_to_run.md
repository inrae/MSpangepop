# How to Run MSpangepop

## Prerequisites

Before running MSpangepop, ensure you have Snakemake installed and access to Singularity containers. The workflow is designed to run on high-performance computing clusters with SLURM scheduling, though it can be adapted for other environments.

## Configuration Setup

Start by editing the master configuration file to define your samples. Each sample requires a compressed FASTA file, population genetics parameters, and simulation settings. The workflow supports multiple samples that can be processed in parallel.

## Resource Considerations

The workflow is computationally intensive, particularly for large genomes or high mutation rates. Memory requirements scale with population size and genome length. The configuration includes a memory multiplier to adjust resource allocation across all workflow steps.

## Execution Environment

The workflow uses Singularity containers to ensure reproducible execution across different computing environments. Container images are pulled from the INRAe registry and contain all required dependencies including msprime, tskit, and visualization libraries.

## Monitoring Progress

The workflow provides detailed logging throughout execution, with color-coded status messages indicating progress through different stages. Each major step reports completion status and timing information to help monitor workflow progress.

## Output Organization

Results are organized in a structured directory hierarchy under the specified output directory. Each sample gets its own subdirectory containing simulation data, graphs, visualizations, and summary files.