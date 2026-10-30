@workflow
Feature: Robustness of the MSpangenome workflow
  These scenarios describe stable behavior and data contracts across the workflow.
  Randomized behavior is checked with fixed seeds or invariant properties, rather
  than relying on exact stochastic output unless a fixed seed is under test.

  @unit @graph
  Scenario Outline: Path operations preserve graph invariants
    Given a deterministic graph fixture with segments "A", "B", and "C"
    And oriented edges forming the path "A-B-C"
    When I apply the path operation "<operation>" to the graph fixture
    Then the resulting path has the expected sequence and orientations
    And every edge in the resulting path references an existing segment

    Examples:
      | operation |
      | bypass    |
      | loop      |
      | invert    |
      | swap      |
      | paste     |
      | cut_paste |

  @unit @graph @validation
  Scenario: Invalid path bounds are rejected clearly
    Given a deterministic graph fixture with a valid path
    When I apply a path operation with a bound outside a segment
    Then the operation raises a clear validation error

  @unit @graph
  Scenario: Reverse complement handles IUPAC bases and letter case
    Given a test sequence containing IUPAC bases and lowercase bases
    When I calculate its reverse complement
    Then the result matches the expected reverse complement
    And non-mutable bases such as "N" follow the documented behavior

  @unit @randomness
  Scenario: A fixed seed makes variant assignment reproducible
    Given the same variant input and random seed "12345"
    When I assign variants twice
    Then both assignments are identical

  @unit @randomness @variants
  Scenario: Variant coordinates and counts respect configured limits
    Given a locus of length "20" and a configured variant limit
    When I assign variants using a fixed seed
    Then every variant position and length is within the locus bounds
    And the number of variants does not exceed the configured limit

  @unit @configuration @randomness
  Scenario: Supported parameter distributions stay within their configured bounds
    Given a valid range for a supported parameter distribution
    When I sample a value using a fixed seed
    Then the sampled value is within the configured bounds

  @unit @configuration
  Scenario Outline: Seed conversion handles supported input forms
    Given the seed value "<seed>"
    When I convert the value to a usable random seed
    Then the conversion result is classified as "<result>"

    Examples:
      | seed       | result                |
      | 12345      | integer               |
      | text-seed  | deterministic_integer |
      | none       | random                |
      | 0          | integer               |

  @unit @configuration @variants
  Scenario: Variant probabilities are normalized and max_sv limits are honored
    Given a valid variant distribution and a max_sv limit for each variant type
    When I normalize the distribution and assign variants using a fixed seed
    Then the normalized probabilities sum to one
    And the count of each variant type does not exceed its max_sv limit

  @integration @traversal
  Scenario: A small tree JSON input produces the expected traversal
    Given a small tree JSON fixture
    When I compute the traversal for the tree
    Then the expected nodes and terminal lineages are present
    And the affected descendants match the fixture

  @integration @variants
  Scenario: An annotated traversal produces variants matching the data contract
    Given a traversal fixture with a few nodes and intervals
    When I generate annotated variants using a fixed seed
    Then the output JSON matches the expected schema
    And variant coordinates, lengths, and limits are valid

  @integration @traversal @graph
  Scenario: A locus containing one node produces a valid graph path
    Given a traversal fixture containing a locus with one node
    When I build the corresponding graph
    Then the locus path contains the expected node
    And every path reference is valid

  @integration @validation
  Scenario: Incompatible lineages or intervals produce a clear error
    Given a fixture with incompatible lineages or interval bounds
    When I pass the fixture to the relevant validation step
    Then a clear error identifies the incompatible data

  @integration @fasta
  Scenario: Intervals extract the expected sequences from a FASTA file
    Given a short FASTA fixture and known intervals
    When I split the FASTA using those intervals
    Then the extracted sequences match the expected coordinates
    And the output headers identify the correct loci

  @integration @gfa
  Scenario: A traversal and its loci produce a consistent GFA graph
    Given an annotated traversal fixture and its corresponding locus sequences
    When I build the graph in GFA format
    Then every segment referenced by a link exists
    And segment identifiers are unique
    And links, orientations, and paths are consistent

  @integration @gfa
  Scenario: Concatenating subgraphs preserves identifiers and path continuity
    Given two small GFA subgraph fixtures with shared lineages
    When I concatenate the subgraphs
    Then segment identifiers are unique in the merged graph
    And inter-subgraph links have the expected orientation
    And paths for shared lineages remain continuous

  @integration @batching
  Scenario Outline: Loci are divided into batches of the configured size
    Given "<locus_count>" loci and a batch size of "<batch_size>"
    When I divide the loci into batches
    Then I get "<batch_count>" batches
    And each batch contains at most "<batch_size>" loci

    Examples:
      | locus_count | batch_size | batch_count |
      | 1           | 100        | 1           |
      | 205         | 100        | 3           |
      | 200         | 100        | 2           |

  @functional @slow @manual_or_scheduled
  Scenario: A small end-to-end simulation produces a usable GFA
    Given a test configuration with few samples and a short chromosome
    And a test FASTA and the required external tools are available
    When I run the workflow through to the final GFA output
    Then the expected output files are present
    And the GFA can be parsed and satisfies its structural invariants
    And every expected lineage has a path
