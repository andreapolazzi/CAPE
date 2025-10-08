# CAPE — Computational Assistant for Pathway Engineering

> A modular bioinformatics tool for pathway prediction, enzyme selection, and gene optimization.

## Overview

**CAPE** is a modular computational platform designed to assist metabolic pathway engineering.
It supports **pathway retrieval** (KEGG) and **prediction** (RetroPath2.0), **EC number suggestion** (SelenzymeRF), **protein sequences fetching** (NCBI Protein), and **codon optimization** for a target host organism.
Each module can be used independently or combined into a full pipeline.

In addition, the whole Django app can be locally hosted with a user-friendly interface.

## Quick Start

### Local host

To run a local host of the web app, you must have both [Django](https://docs.djangoproject.com/en/5.2/topics/install/) and [conda](https://docs.conda.io/en/latest/miniconda.html) installed.

Then, on a Linux-like terminal:

- Clone the GitHub repository

```bash
git clone https://github.com/andreapolazzi/CAPE.git
```
All commands to use this repository are intended to be run from `CAPE/` and directories must be relative to `CAPE/`.

- Create the conda environment (inside the repository)

```bash
cd CAPE
conda env create -f environment.yml
```

- Unzip necessary files
```bash
gunzip -k cape_app/algorithms/pathways/data/retrorules_rr02_rp2_flat_forward.csv.gz
tar -xzf cape_app/algorithms/pathways/data/Pathways.tar.gz -C cape_app/algorithms/pathways/data/
```

- Install RetroPath2.0
```bash
conda install -c conda-forge -n cape retropath2_wrapper
```
Disclaimer: as of 08/10/2025, there is a bug in KNIME download attempt by RetroPath2.0 at its first download. Authors are working to fix the problem, but we cannot guarantee this to work.

- Activate the environment

```bash
conda activate cape
```

- Run the server

```bash
python manage.py runserver
```

You can access the web tool by navigating to <http://127.0.0.1:8000/>.

### Installation

If you want to just use the features of the tool without the web-based implementation, follow the instruction as above (Django is not required in this case).

Each module can be used independently from one another, and they can be streamlined by the user as desired.
See the documentation relative to each module's CLI below.

## Modules

| Step | Description | Link |
|------|--------------|------|
| **1️. Pathway Prediction** | Predicts viable metabolic routes from a starting compound to a desired product. | [Docs →](docs/pathway_prediction.md) |
| **2️. Enzyme Selection** | Retrieves enzyme candidates from databases (MetaNetX, Brenda, EXpasy) to aid identifying missing EC. | [Docs →](docs/enzyme_selection.md) |
| **3️. Protein Retrieval** | Retrieves and exports protein sequences in FASTA format. | [Docs →](docs/protein_retrieval.md) |
| **4️. Codon Optimization** | Adapts coding sequences for the selected chassis organism. | [Docs →](docs/codon_optimization.md) |

As indicated in the documentation files, the scripts can be found under `cape_app/algorithms`.

In the `export/` folder you can find the **HERO parts** and the **_pLoxship_** files, displayed and downloadable in the final export section (for the web tool).

## Pipeline Overview

![Pipeline Diagram](docs/pipeline_diagram.svg)
