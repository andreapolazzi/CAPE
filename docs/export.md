# 🧩 Export Module

Exports optimized sequences for downstream cloning or modeling.

## Supported Formats

- FASTA (default)
- SBOL v3 (via pySBOL3)

## Usage

```bash
python -m cape.export --input optimized.fasta --format sbol
```

## Notes

- Each exported construct includes metadata (gene name, enzyme, EC number).
- SBOL files contain standardized component definitions for interoperability.

[← Back to Home](../README2.md)
