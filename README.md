`test.py` pour créer les arbres correspondants aux chromosomes. Prend un samtools FAI en entrée, donne un fichier trees et un fichier VCF.

`fusion.py` pour fusionner le VCF de msprime, le BED de VISOR et obtenir un VCF final utilisable pour vg.

À faire :
- ajouter l'entête au VCF final (reprendre celui du VCF de msprime)
- donner le nombre de variants dans le VCF msprime --> pour créer un BED avec le même nombre de variants
- ajouter les types de variants manquants dans le script `bed2vcf_2.py`