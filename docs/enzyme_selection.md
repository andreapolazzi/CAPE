# 🧩 Enzyme Selection Module

Identifies enzyme candidates for each predicted reaction.

## Input

- Reaction SMARTS (from pathway prediction)
- Optional EC numbers

## Output

- Table of enzyme candidates (CSV / DataFrame)
- FASTA sequences for selected entries

## Usage

```bash
python -m cape.enzyme_selection --input pathway.csv --output enzymes.csv
```

## Implementation

- Queries [Selenzyme API](http://selenzymerf.synbiochem.co.uk/)
  and falls back to local scoring if unavailable.

## Example Output

| EC | Enzyme | Organism | Score |
|----|---------|-----------|--------|
| 1.1.1.1 | Alcohol dehydrogenase | *E. coli* | 0.91 |

[← Back to Home](../README.md)
