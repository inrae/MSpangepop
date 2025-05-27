# Nucleotide transition matrices for the graph creation

# You can modify and add a custom one

# BLOSUM-INSPIRED MATRIX (for realistic SNP generation)
# Purpose: Generate "realistic" SNPs
# Based on: BLOSUM-like scoring principles adapted for nucleotides
# Key features:
#   - Transition bias: A<->G and C<->T mutations are more common
#   - Transversion penalty: A<->C, A<->T, G<->C, G<->T mutations are less common
blosum_matrix = {
    "A": {"C": 0.25, "G": 0.50, "T": 0.25},  # A->G transition favored (0.50)
    "C": {"A": 0.30, "G": 0.20, "T": 0.50},  # C->T transition favored (0.50)  
    "G": {"A": 0.50, "C": 0.20, "T": 0.30},  # G->A transition favored (0.50)
    "T": {"A": 0.25, "C": 0.50, "G": 0.25}   # T->C transition favored (0.50)
}


# RANDOM MATRIX (for neutral SNP generation)
# Purpose: Generate completely random single nucleotide polymorphisms (SNPs)
# Based on: Jukes-Cantor model assumption of equal mutation rates
# Key features:
#   - No mutational bias: all substitutions equally likely (1/3 probability each)
random_matrix = {
    "A": {"C": 1/3, "G": 1/3, "T": 1/3},
    "C": {"A": 1/3, "G": 1/3, "T": 1/3},
    "G": {"A": 1/3, "C": 1/3, "T": 1/3},  
    "T": {"A": 1/3, "C": 1/3, "G": 1/3}   
}


# AT-BIASED MATRIX (for transposable element sequence generation)
# Purpose: Generate AT-rich sequences mimicking transposable elements (TEs)
# Based on: Empirical observations that TEs are typically 60-70% AT content (HKY85)
# Key features:
#   - AT bias: A and T are more likely to transition to each other than to G/C
#   - Reflects TE compositional bias found in many organisms
simple_at_bias_matrix = {
    '_' : {'A': 0.70, 'T': 0.20, 'G': 0.06, 'C': 0.04}, # Etat initial 
    'A': {'A': 0.70, 'T': 0.20, 'G': 0.06, 'C': 0.04},  # Strong A retention, bias toward T
    'T': {'A': 0.20, 'T': 0.70, 'G': 0.04, 'C': 0.06},  # Strong T retention, bias toward A  
    'G': {'A': 0.25, 'T': 0.15, 'G': 0.45, 'C': 0.15},  # Moderate G retention, AT bias
    'C': {'A': 0.15, 'T': 0.25, 'G': 0.15, 'C': 0.45}   # Moderate C retention, AT bias
}

# TEST MATRIX
# Purpose: Be able to locate insertions in the graph
test_matrix = {
    '_' : {'A': 1, 'T': 0, 'G': 0, 'C': 0},
    'A': {'A': 1, 'T': 0, 'G': 0, 'C': 0},
    'T': {'A': 1, 'T': 0, 'G': 0, 'C': 0},
    'G': {'A': 1, 'T': 0, 'G': 0, 'C': 0},
    'C': {'A': 1, 'T': 0, 'G': 0, 'C': 0}
}