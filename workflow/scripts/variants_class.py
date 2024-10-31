import random

class Variant:
    """Base class for genetic variants."""
    bases = ['A', 'C', 'G', 'T']
    repetition_range = (2, 3)

    def __init__(self, variant_type, chrom, position, chromosome_dict, chrom_lengths, length=1, reference_seq="", samples=None):
        self.variant_type = variant_type
        self.chrom = chrom
        self.position = position
        self.chromosome_dict = chromosome_dict
        self.chrom_lengths = chrom_lengths
        self.length = length
        self.reference_seq = reference_seq
        self.alt_seq = ""
        self.samples = samples if samples is not None else {}

    def set_length(self, length_df):
        """Determine variant length based on provided length distribution."""
        rand_val = random.random()
        row = length_df[length_df['cumulative_pb'] >= rand_val].iloc[0]
        interval = row['size_interval']
        lower_bound, upper_bound = map(float, interval.strip('[]').split(','))
        self.length = random.randint(int(lower_bound), int(upper_bound))

        # Ensure length does not exceed chromosome boundaries
        max_length = self.chrom_lengths[self.chrom] - self.position
        self.length = min(self.length, max_length)
        if self.length == 0:
            self.length = 1

    def retrieve_reference_sequence(self, chrom=None, position=None, no_ref=False):
        """Retrieve the reference sequence from a FASTA file."""
        chrom = chrom or self.chrom
        position = position or self.position
        
        if chrom in self.chromosome_dict:
            seq_record = self.chromosome_dict[chrom]
            start = position - 1
            if no_ref:
                return str(seq_record.seq[start:start + 1])  # Only one base needed for insertions
            end = min(start + self.length, len(seq_record))
            return str(seq_record.seq[start:end])
        else:
            raise ValueError(f"Chromosome {chrom} not found in FASTA.")

    def compute_alt_seq(self):
        raise NotImplementedError("Subclasses must implement this method.")

    def reverse_sequence(self, sequence):
        return sequence[::-1]

    def describe(self):
        return f"{self.variant_type}\t{self.chrom}\t{self.position}\t{self.length}\t{self.reference_seq[:10]}\t{self.alt_seq[:10]}\t{self.samples}"

    def vcf_line(self):
        samples_str = '\t'.join(map(str, self.samples.values()))
        return f"{self.chrom}\t{self.position}\t.\t{self.reference_seq}\t{self.alt_seq}\t.\tPASS\t.\tGT\t{samples_str}"

class SingleBaseVariant(Variant):
    """Variants that involve a change at a single base."""
    
    def __init__(self, variant_type, chrom, position, chromosome_dict, chrom_lengths):
        super().__init__(variant_type, chrom, position, chromosome_dict, chrom_lengths, length=1)

class SNP(SingleBaseVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths):
        super().__init__('SNP', chrom, position, chromosome_dict, chrom_lengths)

    def compute_alt_seq(self):
        self.alt_seq = random.choice([base for base in self.bases if base != self.reference_seq])

class StructuralVariant(Variant):
    """Variants that involve multiple bases and structural changes."""
    
    def __init__(self, variant_type, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__(variant_type, chrom, position, chromosome_dict, chrom_lengths)
        self.set_length(length_df)

class Deletion(StructuralVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('Deletion', chrom, position, chromosome_dict, chrom_lengths, length_df)

    def compute_alt_seq(self):
        self.alt_seq = self.reference_seq[0]  # Simulate deletion

class Insertion(StructuralVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('Insertion', chrom, position, chromosome_dict, chrom_lengths, length_df)
        self.retrieve_reference_sequence(no_ref=True)

    def compute_alt_seq(self):
        self.alt_seq = self.reference_seq + "".join(random.choice(self.bases) for _ in range(self.length))

class Inversion(StructuralVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('Inversion', chrom, position, chromosome_dict, chrom_lengths, length_df)

    def compute_alt_seq(self):
        self.alt_seq = self.reverse_sequence(self.reference_seq)

class RepeatableVariant(StructuralVariant):
    """Variants that involve repetitions."""
    
    def compute_alt_seq(self):
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reference_seq * number_of_repeats

class TandemDuplication(RepeatableVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('TandemDuplication', chrom, position, chromosome_dict, chrom_lengths, length_df)

class InvertedTandemDuplication(RepeatableVariant):
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('InvertedTandemDuplication', chrom, position, chromosome_dict, chrom_lengths, length_df)

    def compute_alt_seq(self):
        number_of_repeats = random.randint(self.repetition_range[0], self.repetition_range[1])
        self.alt_seq = self.reverse_sequence(self.reference_seq * number_of_repeats)

class Translocation(Variant):
    """Variants that involve translocations."""
    def __init__(self, chrom, position, chromosome_dict, chrom_lengths, length_df):
        super().__init__('Translocation', chrom, position, chromosome_dict, chrom_lengths)
        self.set_length(length_df)
        self.dest_chrom = None
        self.dest_position = None
        self.dest_reference_seq = ""  # This must be one base at dest_position on dest_chrom
        self.dest_alt_seq = ""

    def set_destination(self, dest_chrom, dest_position):
        """Set the destination chromosome and position, retrieving the reference sequence."""
        self.dest_chrom = dest_chrom
        self.dest_position = dest_position
        
        # Retrieve the reference sequence for the destination position
        self.dest_reference_seq = self.retrieve_reference_sequence(self.dest_chrom, self.dest_position, no_ref=True)
        self.dest_alt_seq = ""  # Initialize the destination alternative sequence

    def compute_alt_seq(self):
        """Compute the alternative sequences for the translocation."""
        if self.dest_chrom is not None and self.dest_position is not None:
            # Set the alternative sequence for the original position
            self.alt_seq = self.reference_seq[0]  # Use the first base from the original sequence
            # Set the destination alternative sequence to be the same as the destination reference sequence
            self.dest_alt_seq = self.dest_reference_seq + self.reference_seq

    def describe(self):
        """Return a description of the translocation variant, including destination info."""
        return f"{self.variant_type}\t{self.chrom}\t{self.position}\t{self.dest_chrom}\t{self.dest_position}\t{self.dest_alt_seq[:10]}\t{self.samples}"

    def vcf_line(self):
        """Generate a VCF line for both the original and destination variants."""
        samples_str = '\t'.join(map(str, self.samples.values()))
        orig_vcf = f"{self.chrom}\t{self.position}\t.\t{self.reference_seq}\t{self.alt_seq}\t.\tPASS\t.\tGT\t{samples_str}"
        dest_vcf = f"{self.dest_chrom}\t{self.dest_position}\t.\t{self.dest_reference_seq}\t{self.dest_alt_seq}\t.\tPASS\t.\tGT\t{samples_str}"
        return f"{orig_vcf}\n{dest_vcf}"
