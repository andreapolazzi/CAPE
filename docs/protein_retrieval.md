# 🧩 Protein Retrieval Module

Fetches protein sequences corresponding to enzyme candidates.

## Input

- EC numbers or accession IDs from enzyme selection output.

## Output

- FASTA file with protein sequences.

## Usage

```bash
python -m cape.protein_retrieval --input enzymes.csv --output proteins.fasta
```

## Implementation

- Uses NCBI Entrez API for sequence retrieval.
- Returns all matching entries for user curation.

[← Back to Home](../README2.md)
