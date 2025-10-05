# 🧩 Codon Optimization Module

Optimizes coding sequences for expression in the target chassis organism.

## Input

- FASTA file of protein sequences
- Host organism name or codon usage table

## Output

- Codon-optimized DNA sequences (FASTA)

## Usage

```bash
python -m cape.codon_optimization --input proteins.fasta --host rhodococcus --output optimized.fasta
```

## Implementation

- Respects host codon usage bias.
- Avoids restriction sites (BioBrick / Golden Gate compatible).

[← Back to Home](../README2.md)
