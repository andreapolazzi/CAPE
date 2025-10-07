# Pathway Prediction Module

The pathway prediction logic resides in:
`cape_app/algorithms/pathways/main.py`

Helper classes:

- `cape_app/algorithms/pathways/classes/chem_utils.py`
- `cape_app/algorithms/pathways/classes/classes_functions.py`
- `cape_app/algorithms/pathways/classes/rp_classes.py`

This module identifies plausible metabolic routes between a **source** and a **target compound**, integrating both curated KEGG reactions and *de novo* routes generated via **RetroPath2.0**.  
It constructs a weighted metabolic graph, computes shortest or k-shortest paths, and dynamically extends the network when no direct route is found.

The module can be invoked from the command line through its built-in CLI.

---

## Input

- Starting and target compound InChI strings (mandatory).
- Optional compound names for readable labels in output.
- Optional RetroPath2.0 configuration (reaction rule files, sink file, etc.).

---

## Output

- CSV file (`pathways.csv`) containing the k-shortest paths found between the input compounds.
- Each row lists compound identifiers, reaction identifiers, and path metadata.
- Output directory defaults to `cape_app/pipeline_results/rp_results/`.

---

## Usage

```bash
# Run a standard pathway search using KEGG + RetroPath integration
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"   --end_inchi "InChI=1S/C2H4O/c1-2-3/h2H,1H3"

# Specify a custom output folder and save the results as CSV
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C7H8/c1-2-4-6-7-5-3-1/h1-7H"   --end_inchi "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"   --outdir results/pathway_search

# Tune RetroPath parameters for deeper exploration
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C8H10O/c1-9-8-6-4-2-3-5-7-8/h2-7H,1H3"   --end_inchi "InChI=1S/C4H6O2/c1-3(5)4(2)6/h1H,2H3"   --max_paths 5   --max_steps 6   --max_structures 150   --min_diameter 8   --max_diameter 18
```

---

## Implementation

- Builds a **graph of metabolites** using curated KEGG pathways via `build_graph_from_tsv()`.
- Standardises input InChIs to InChIKeys for unambiguous node representation (`Std.standardize`).
- Computes k-shortest paths via Dijkstra-like traversal (`k_shortest_paths`).
- If the goal compound is absent from the network, calls **RetroPath2.0** through `RP.django_call()` to generate new reactions.
- Automatically merges RetroPath edges into the graph and retries path search.
- Produces a `.csv` summary of all retrieved paths with associated weights and compound identifiers.

### Key Arguments

- `--start_inchi`, `--end_inchi`: required; define the endpoints of the metabolic search.
- `--max_paths`: maximum number of shortest paths to return (default 3).
- `--max_steps`: maximum RetroPath iteration depth (default 5).
- `--max_structures`: maximum number of structures generated per step.
- `--min_diameter`, `--max_diameter`: control reaction rule diameter range (default 6–16).
- `--timeout`: execution limit for RetroPath subprocess.
- `--sink_file`, `--rules_file`: optional custom paths to RetroPath input data.
- `--outdir`: output directory; created automatically if missing.

---

## Example Output

| Path ID | Start InChIKey | End InChIKey | Steps | Method     |
|----------|----------------|--------------|--------|-------------|
| 1        | UHOVQNZJYSORNB-UHFFFAOYSA-N | LFQSCWFLJHTTHZ-UHFFFAOYSA-N | 4 | KEGG + RetroPath |
| 2        | UHOVQNZJYSORNB-UHFFFAOYSA-N | LFQSCWFLJHTTHZ-UHFFFAOYSA-N | 5 | RetroPath only   |

Saved to:  
`cape_app/pipeline_results/rp_results/pathways.csv`

---

## Customisation Tips

- **Different rule sets**: swap `--rules_file` with any RetroRules export to target specific reaction classes.
- **Partial KEGG networks**: edit the local TSV files to restrict search space to certain organisms or pathways.
- **Heavier exploration**: increase `--max_structures` and `--max_steps` for complex transformations (at the cost of runtime).
- **Batch usage**: combine multiple runs in shell loops or workflow managers (e.g. Snakemake).

---

[← Back to Home](../README.md)
