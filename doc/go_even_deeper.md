# Going Even Deeper

## Readable JSON

MSpangepop includes an option to output JSON files in human-readable format with proper indentation and formatting. While this makes the files larger and slower to process, it enables manual inspection and debugging.

## Markov Matrices for Insertions

The workflow uses Markov chain models to generate realistic DNA sequences for insertions. These transition matrices can be customized to produce sequences with specific compositional biases or evolutionary characteristics. You can find them in the `workflow/scripts/matrix.py`

## Custom Mutation Models

The workflow supports multiple mutation models beyond the default binary model, including :
- infinite alleles
- Jukes-Cantor models

## Scalability Considerations

MSpangepop is designed to scale from small test datasets to whole-genome simulations with thousands of individuals. Resource requirements and processing strategies automatically adapt to dataset size and available computational resources.
