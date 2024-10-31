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
        """Must be defined in subclasses to compute the alternative sequence."""
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


class SingleBaseVariant(Variant):
    """Intermediate class for single-base modifications."""
    def __init__(self, variant_type, chrom, position):
        super().__init__(variant_type, chrom, position, length=1)


class SNP(SingleBaseVariant):
    def __init__(self, chrom, position):
        super().__init__('SNP', chrom, position)

    def compute_alt_seq(self):
        """Creates an alternative sequence by changing a single base."""
        self.alt_seq = random.choice([base for base in self.bases if base != self.reference_seq])


class StructuralVariant(Variant):
    """Intermediate class for variants that involve length-based sequence modifications."""
    def __init__(self, variant_type, chrom, position, length):
        super().__init__(variant_type, chrom, position, length)


class Deletion(StructuralVariant):
    def __init__(self, chrom, position, length):
        super().__init__('Deletion', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by simulating a deletion."""
        self.alt_seq = self.reference_seq[0]  # Assumes deletion keeps one base


class Insertion(StructuralVariant):
    def __init__(self, chrom, position, length):
        super().__init__('Insertion', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by adding random bases."""
        self.alt_seq = self.reference_seq + "".join(random.choice(self.bases) for _ in range(self.length))


class Inversion(StructuralVariant):
    def __init__(self, chrom, position, length):
        super().__init__('Inversion', chrom, position, length)
    
    def compute_alt_seq(self):
        """Creates an alternative sequence by reversing the reference sequence."""
        self.alt_seq = self.reverse_sequence(self.reference_seq)


class RepeatableVariant(Variant):
    """Intermediate class for variants that involve repeated sequences."""
    def __init__(self, variant_type, chrom, position, length):
        super().__init__(variant_type, chrom, position, length)

    def compute_repeated_sequence(self, reference_seq):
        """Creates a repeated sequence based on the repetition range."""
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        return reference_seq * number_of_repeats


class TandemDuplication(RepeatableVariant):
    def __init__(self, chrom, position, length):
        super().__init__('TandemDuplication', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by duplicating the reference sequence."""
        self.alt_seq = self.compute_repeated_sequence(self.reference_seq)


class InvertedTandemDuplication(RepeatableVariant):
    def __init__(self, chrom, position, length):
        super().__init__('InvertedTandemDuplication', chrom, position, length)

    def compute_alt_seq(self):
        """Creates an alternative sequence by duplicating and reversing the reference sequence."""
        self.alt_seq = self.reverse_sequence(self.compute_repeated_sequence(self.reference_seq))


class Translocation(StructuralVariant):
    def __init__(self, chrom, position, length):
        super().__init__('Translocation', chrom, position, length)

    def compute_alt_seq(self):
        """Placeholder for translocation logic, if applicable."""
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