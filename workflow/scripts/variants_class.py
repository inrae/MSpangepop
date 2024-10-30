"""
Author: Lucien Piat
Date: 28 Oct 2024
Institution: INRAe
Project: PangenOak
"""

import random

class Variant:
    """Base class for genetic variants."""
    bases = ['A', 'C', 'G', 'T']
    repetition_range = (2, 3)

    def __init__(self, variant_type, chrom, pos, length=1, reference_seq="", samples=None):
        self.variant_type = variant_type
        self.chrom = chrom
        self.pos = pos
        self.length = length
        self.reference_seq = reference_seq
        self.alt_seq = ""
        self.samples = samples if samples is not None else {}

    def compute_alt_seq(self):
        #Must be defined in subclasses
        pass 

    def reverse_sequence(self, sequence):
        return sequence[::-1]

    #This is the line in the recap file
    def describe(self):
        return f"{self.variant_type}\t{self.chrom}\t{self.pos}\t{self.length}\t{self.reference_seq[:10]}\t{self.alt_seq[:10]}\t{self.samples}"

    #This is the line printed in the vcf
    def vcf_line(self):
        samples_str = '\t'.join(map(str, self.samples.values()))
        return f"{self.chrom}\t{self.pos}\t.\t{self.reference_seq}\t{self.alt_seq}\t.\tPASS\t.\tGT\t{samples_str}"

class SNP(Variant):
    def __init__(self, chrom, pos):
        super().__init__('SNP', chrom, pos, length=1)

    def compute_alt_seq(self):
        self.alt_seq = random.choice([base for base in self.bases if base != self.reference_seq])
        
class Deletion(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('Deletion', chrom, pos, length)

    def compute_alt_seq(self):
        self.alt_seq = self.reference_seq[0]

class Insertion(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('Insertion', chrom, pos, length)

    def compute_alt_seq(self):
        self.alt_seq = self.reference_seq + "".join(random.choice(self.bases) for _ in range(self.length))

class Inversion(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('Inversion', chrom, pos, length)
    
    def compute_alt_seq(self):
        self.alt_seq = self.reverse_sequence(self.reference_seq)

class TandemDuplication(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('TandemDuplication', chrom, pos, length)

    def compute_alt_seq(self):
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reference_seq * number_of_repeats
        
class InvertedTandemDuplication(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('InvertedTandemDuplication', chrom, pos, length)
    
    def compute_alt_seq(self):
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reverse_sequence(self.reference_seq * number_of_repeats)

class Translocation(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('Translocation', chrom, pos, length)

class Transduplication(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('Transduplication', chrom, pos, length)

class ReciprocalTranslocation(Variant):
    def __init__(self, chrom, pos, length):
        super().__init__('ReciprocalTranslocation', chrom, pos, length)