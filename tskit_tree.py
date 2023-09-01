import tskit

ts = tskit.load("combined.trees")  # Or generate using e.g. msprime.sim_ancestry()
print(ts)  # In a Jupyter notebook this displays a summary table. Otherwise use print(ts)
genetic_diversity = ts.diversity()
print("Av. genetic diversity across the genome is", genetic_diversity)

t = ts.as_vcf()
print(t)