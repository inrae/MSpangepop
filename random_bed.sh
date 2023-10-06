### create random BED

# refFAI="ztIPO323.chr8.fasta.fai"

# cut -f1,2 $refFAI > chrom.dim.tsv
# Rscript randomregion.r -d chrom.dim.tsv -n 2638 -v 'deletion,inversion,inverted tandem duplication,translocation copy-paste,reciprocal translocation' -r '40:20:20:10:10' | sortBed > HACk.random.bed


### vg construct
IMG="/home/sukanya/singularity/builddir/vg_v1.51.0.sif"
REF="/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/ztIPO323.chr8.fasta"
VCF="/home/sukanya/repos/00_to_forgemia/msprime_VISOR_fusionVCF/results/60_tritici.vcf"

# $IMG vg construct -r $REF -v $VCF -S -f -m 250000000 -a >test/test.vg
$IMG vg construct -r $REF -v $VCF -S -f -m 250000000 >test/test.vg

# # $IMG vg convert
# $IMG vg paths -F -v test/test.vg >test/test.fasta
# $IMG vg paths -E -x test/test.vg

$IMG vg convert -f test/test.vg > test/test.gfa

# $IMG vg deconstruct -p chr_8 test/test.gfa > test/test.vcf