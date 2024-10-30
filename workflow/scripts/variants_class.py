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

    def __init__(self, variant_type, chrom, position, length=1, reference_seq="", samples=None):
        self.variant_type = variant_type
        self.chrom = chrom
        self.position = position
        self.length = length
        self.reference_seq = reference_seq
        self.alt_seq = ""
        self.samples = samples if samples is not None else {}

    def compute_alt_seq(self):
        #Must be defined in subclasses
        raise NotImplementedError("Subclasses must implement this method.")

    def reverse_sequence(self, sequence):
        """Reverses a given sequence."""
        return sequence[::-1]

    def describe(self):
        """Creates a summary line for this variant."""
        return f"{self.variant_type}\t{self.chrom}\t{self.position}\t{self.length}\t{self.reference_seq[:10]}\t{self.alt_seq[:10]}\t{self.samples}"

    def vcf_line(self):
        """Formats the variant as a line in VCF format."""
        samples_str = '\t'.join(map(str, self.samples.values()))
        return f"{self.chrom}\t{self.position}\t.\t{self.reference_seq}\t{self.alt_seq}\t.\tPASS\t.\tGT\t{samples_str}"

class SNP(Variant):
    def __init__(self, chrom, position):
        super().__init__('SNP', chrom, position, length=1)

    def compute_alt_seq(self):
        """Creates an alternative sequence by changing a single base."""
        self.alt_seq = random.choice([base for base in self.bases if base != self.reference_seq])
        
class Deletion(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('Deletion', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by keeping only the first base, simulating a deletion."""
        self.alt_seq = self.reference_seq[0]

class Insertion(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('Insertion', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by adding random bases."""
        self.alt_seq = self.reference_seq + "".join(random.choice(self.bases) for _ in range(self.length))

class Inversion(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('Inversion', chrom, position, length)
    
    def compute_alt_seq(self):
        """Creates an alternative sequence by reversing the reference sequence."""
        self.alt_seq = self.reverse_sequence(self.reference_seq)

class TandemDuplication(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('TandemDuplication', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by duplicating the reference sequence."""
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reference_seq * number_of_repeats
        
class InvertedTandemDuplication(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('InvertedTandemDuplication', chrom, position, length)
    
    def compute_alt_seq(self):
        """Creates an alternative sequence by duplicating and reversing the reference sequence."""
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reverse_sequence(self.reference_seq * number_of_repeats)

class Translocation(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('Translocation', chrom, position, length)
    
    def compute_alt_seq(self):
        """Placeholder method for translocation type variant."""
        pass 

class Transduplication(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('Transduplication', chrom, position, length)
    
    def compute_alt_seq(self):
        """Placeholder method for transduplication type variant."""
        pass 

class ReciprocalTranslocation(Variant):
    def __init__(self, chrom, position, length):
        super().__init__('ReciprocalTranslocation', chrom, position, length)

    def compute_alt_seq(self):
        """Placeholder method for reciprocal translocation type variant."""
        pass