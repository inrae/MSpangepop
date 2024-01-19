
### outside of msprime_VISOR_fusionVCF folder
sudo singularity build mspvisor.sif msprime_VISOR_fusionVCF/singularity_img/to_singularity.def

### use singularity image
# singularity exec mspvisor.sif python3 /mspv/tree_generation.py --help
