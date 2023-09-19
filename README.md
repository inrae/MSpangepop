Depuis le script `main.py`:
- génération d'une population avec msprime : `msprime_vcf()`
    - input : FASTA index, taille de pop, mutation rate, nombre d'individus à simuler
    - output : `msprime_output_chr.vcf`
- fusion du BED de VISOR et du VCF de msprime : `main()`
    - input : BED pour VISOR, VCF msprime, FASTA, nom du fichier output
    - output : VCF

`tree_generation.py` pour créer les arbres correspondants aux chromosomes. Prend un samtools FAI en entrée, donne un fichier trees et un fichier VCF.

`fusion.py` pour fusionner le VCF de msprime, le BED de VISOR et obtenir un VCF final utilisable pour vg.

À faire :
- ajouter l'entête au VCF final (reprendre celui du VCF de msprime)
- donner le nombre de variants dans le VCF msprime --> pour créer un BED avec le même nombre de variants
- ajouter les types de variants manquants dans le script `bed2vcf.py`