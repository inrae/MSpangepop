import random
import itertools
from fasta_to_gfa import Node
from matrix import random_matrix as transition_matrix

class Variant:
    """Base class for genetic variants."""
    bases = "ATCG"

    def __init__(self, A: Node, B: Node = None, C: Node = None, D: Node = None, length: int = 0, transition_matrix: dict = transition_matrix):
        self.length: int = length
        self.A: Node = A  # Previous node
        self.B: Node = B  # Current node (affected node)
        self.C: Node = C  # Another affected node (if needed)
        self.D: Node = D  # End node
        self.transition_matrix = transition_matrix

    def compute_alt_seq(self, node_id_generator: itertools.count):
        raise NotImplementedError("Subclasses must implement this method.")

class SNP(Variant):
    def __init__(self, A: Node, B: Node, D: Node, transition_matrix: dict = transition_matrix):
        """
        SNP modifies the nucleotide of the `B` node based on the transition matrix.
        """
        super().__init__(A, B, None, D, 1, transition_matrix)

    def compute_alt_seq(self, node_id_generator: itertools.count) -> Node:
        """Creates an alternative nucleotide at the SNP position based on the transition matrix."""
        original_base = self.B.base.decode()
        new_base = self._mutate_base(original_base)

        # Create new SNP node
        new_node = Node(new_base, next(node_id_generator))

        # Update graph structure
        self.A.connect(new_node)
        new_node.connect(self.D)

        return new_node

    def _mutate_base(self, original_base: str) -> str:
        """Uses the provided transition matrix to determine the mutated base."""
        return random.choices(
            population=list(self.transition_matrix[original_base].keys()),
            weights=list(self.transition_matrix[original_base].values()),
            k=1
        )[0]

class Deletion(Variant):
    def __init__(self, A: Node, D: Node):
        super().__init__(A, None, None, D, 0)

    def compute_alt_seq(self) -> None:
        self.A.connect(self.D)  # Bypass intermediate nodes

class Insertion(Variant):
    def __init__(self, A: Node, D: Node, length: int):
        super().__init__(A, None, None, D, length)

    def compute_alt_seq(self, node_id_generator: itertools.count) -> list:
        """Generates a list of random one-base nodes and inserts them between A and D."""
        inserted_nodes = []

        # Generate 'length' number of random base nodes
        for _ in range(self.length):
            new_base = random.choice(self.bases)  # Pick a random nucleotide
            new_node = Node(new_base, next(node_id_generator))
            inserted_nodes.append(new_node)

        # Update graph connections
        self.A.connect(inserted_nodes[0])

        for i in range(len(inserted_nodes) - 1):
            inserted_nodes[i].connect(inserted_nodes[i + 1])

        inserted_nodes[-1].connect(self.D)

        return inserted_nodes

class Inversion(Variant):
    """Represents an inversion mutation in the graph."""
    
    def __init__(self, A: Node, B: Node, C: Node, D: Node):
        """Initializes an Inversion variant, flipping the order of B and C."""
        super().__init__(A, B, C, D, 0)

    def compute_alt_seq(self) -> None:
        """Reorders nodes B and C in reverse orientation while maintaining connections."""
        self.A.connect(self.C)
        self.C.connect(self.B)
        self.B.connect(self.C)
        self.B.connect(self.D)
    # TODO change this logic