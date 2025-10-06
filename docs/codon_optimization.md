# Codon Optimization Module

The script is located at:
`cape_app/algorithms/codon_opt/backtranslate_pd630.py`

The codon usage table used can be found at <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37919>

This module optimizes coding sequences for expression in the target chassis organism and prevents the formation of illegal sites (e.g. restriction sites) decided by the user.

## Input

- .faa file of protein sequences

## Output

- Codon-optimized DNA sequences (.fna) using the same FASTA headers as the protein input.

## Usage

```bash
# Deterministic (highest-frequency codons) back-translation
python cape_app/algorithms/codon_opt/backtranslate_pd630.py \
  --input path/to/proteins.faa \
  --output path/to/optimized_sequences.fna \
  --mode max \
  --assembly typeIIS

# Weighted (frequency-proportional) option with a reproducible RNG seed
python cape_app/algorithms/codon_opt/backtranslate_pd630.py \
  --input path/to/proteins.faa \
  --output path/to/optimized_sequences.fna \
  --mode weighted \
  --seed 42 \
  --assembly typeIIS biobrick \
  --forbid CCTGCAGG # additional custom site(s)
```

## Implementation

- Parses protein FASTA files and preserves the original headers in the DNA FASTA output.
- Respects host codon usage bias by ranking synonymous codons using PD630 frequencies (``--mode max`` picks the highest-frequency codon, ``--mode weighted`` samples codons proportional to frequency using ``--seed``).
- Applies any combination of assembly presets (``--assembly``) and user-provided motifs (``--forbid``) to forbid both sense and reverse-complement restriction sites.
- Uses bounded depth-first backtracking (``--max-backtrack``) to resolve conflicts when all codons for a residue would otherwise form illegal sites.
- Formats the resulting DNA with 70-character FASTA line wrapping and reports progress to stdout.

### Key Arguments

- ``--input``/``--output``: required paths to the .faa source and .fna sink.
- ``--mode``: ``max`` (default) or ``weighted`` stochastic sampling; pairing with ``--seed`` ensures reproducibility.
- ``--assembly``: optional preset(s) for common cloning standards (``typeIIS``, ``biobrick``); defaults to ``none``.
- ``--forbid``: extra motifs to avoid beyond the presets; accepts multiple space-separated sequences.
- ``--max-backtrack``: limits how many prior residues can be revisited when searching for a legal codon.

### Customisation Tips

- **Different host organism**: duplicate the script and replace ``CODON_FREQ_PD630`` with host-specific frequencies (Kazusa or literature) while keeping the rest of the pipeline identical.
- **Alternative genetic code**: update ``CODON_TO_AA`` (and regenerate ``AA_TO_CODONS``) if the host uses a non-standard code.
- **Preset list maintenance**: extend ``ASSEMBLY_PRESETS`` for additional cloning schemes; user-level overrides continue to work via ``--assembly``.

[← Back to Home](../README.md)
