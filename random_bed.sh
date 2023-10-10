### create random BED

refFAI="ztIPO323.chr8.fasta.fai"

cut -f1,2 $refFAI > chrom.dim.tsv
Rscript randomregion.r -d chrom.dim.tsv -n 2638 -v 'deletion,inversion,inverted tandem duplication,translocation copy-paste,reciprocal translocation' -r '40:20:20:10:10' | sortBed > HACk.random.bed