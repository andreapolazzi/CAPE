# 🚀 Quickstart Guide

## 1. Install CAPE

```bash
git clone https://github.com/andreapolazzi/CAPE.git
cd CAPE
pip install -r requirements.txt
```

## 2. Prepare Your Input

Provide a `.smi` file containing the substrate, or specify it directly:

```bash
python cape.py --input substrate.smi --host rhodococcus
```

## 3. Run the Full Pipeline

CAPE will:

1. Predict a metabolic pathway.
2. Retrieve enzyme candidates.
3. Fetch corresponding protein sequences.
4. Perform codon optimization.
5. Export results.

Output files (FASTA, CSV, TSV) from each step will be stored in the corresponding folders under `cape_app/`.

See individual module docs for advanced options.

[← Back to Home](../README2.md)
