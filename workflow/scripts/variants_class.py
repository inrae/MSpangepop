import random
import itertools
from fasta_to_gfa import Node
from matrix import random_matrix as transition_matrix

class Variant:
    """Base class for genetic variants."""
    
    bases = "ATCG"

    def __init__(self, start_node: Node, length: int, transition_matrix: dict = transition_matrix):
        """
        Initializes a genetic variant.
        """
        self.length: int = length
        self.start_node: Node = start_node
        self.transition_matrix = transition_matrix

    def compute_alt_seq(self, node_id_generator: itertools.count) -> Node:
        """To be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement this method.")

class SNP(Variant):
    """Represents a Single Nucleotide Polymorphism (SNP)."""
    
    def __init__(self, start_node: Node, transition_matrix: dict = transition_matrix):
        super().__init__(start_node, 1, transition_matrix)

    def compute_alt_seq(self, node_id_generator: itertools.count) -> Node:
        """
        Creates an alternative nucleotide at the SNP position based on the transition matrix.
        """
        original_base = self.start_node.base.decode()
        new_base = self._mutate_base(original_base)

        new_node = Node(new_base, next(node_id_generator))

        # Copy edges from the original node
        new_node.in_edges = list(self.start_node.in_edges)
        new_node.out_edges = list(self.start_node.out_edges)

        # Update connections in the graph
        for in_node in new_node.in_edges:
            in_node.out_edges.append(new_node)
        for out_node in new_node.out_edges:
            out_node.in_edges.append(new_node)

        return new_node

    def _mutate_base(self, original_base: str) -> str:
        """
        Uses the provided transition matrix to determine the mutated base.
        """
        if original_base not in self.transition_matrix:
            raise ValueError(f"Base {original_base} not found in transition matrix.")

        return random.choices(
            population=list(self.transition_matrix[original_base].keys()),
            weights=list(self.transition_matrix[original_base].values()),
            k=1
        )[0]
