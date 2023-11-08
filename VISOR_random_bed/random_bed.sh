### create random BED
refFAI="/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/ztIPO323.chr8.fasta.fai"
GENOME="/home/sukanya/tests/02_data/hackathon_Ztritici/CHR8/ztIPO323.chr8.fasta"
VISOR_script="randomregion.r"

# samtools faidx $GENOME
cut -f1,2 $refFAI > chrom.dim.tsv
Rscript $VISOR_script -d chrom.dim.tsv -n 10 -l 2 -v 'SNP,insertion,deletion,inversion' -r '70:10:10:10' -g $GENOME | sortBed > HACk.random.bed