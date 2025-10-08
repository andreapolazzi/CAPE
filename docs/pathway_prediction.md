# Pathway Prediction Module

The pathway prediction logic resides in:
`cape_app/algorithms/pathways/main.py`

Helper classes:

- `cape_app/algorithms/pathways/classes/chem_utils.py`
- `cape_app/algorithms/pathways/classes/graph_functions.py`
- `cape_app/algorithms/pathways/classes/rp_classes.py`

This module identifies plausible metabolic routes between a **source** and a **target compound**, integrating both curated KEGG reactions and *de novo* routes generated via **RetroPath2.0**.  
It constructs a weighted metabolic graph, computes shortest or k-shortest paths, and dynamically extends the network when no direct route is found.

The module can be invoked from the command line through its built-in CLI.

---

## Input

- Starting and target compound InChI strings (mandatory).
- Optional compound names for readable labels in output.
- Optional number of paths and RetroPath2.0 configuration (reaction rule files, sink file, etc.).

---

## Output

- CSV file (`pathways.csv`) containing the k-shortest paths found between the input compounds.
- Each row represents a reaction, listing compound identifiers, reaction identifiers, and path metadata.
- Output directory defaults to `cape_app/pipeline_results/rp_results/`.

---

## Usage

```bash
# Run a standard pathway search using KEGG + RetroPath integration
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"   --end_inchi "InChI=1S/C2H4O/c1-2-3/h2H,1H3"

# Specify a custom output folder and save the results as CSV
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C7H8/c1-2-4-6-7-5-3-1/h1-7H"   --end_inchi "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"   --outdir results/pathway_search

# Tune RetroPath parameters for deeper exploration
python -m cape_app.algorithms.pathways.main   --start_inchi "InChI=1S/C8H10O/c1-9-8-6-4-2-3-5-7-8/h2-7H,1H3"   --end_inchi "InChI=1S/C4H6O2/c1-3(5)4(2)6/h1H,2H3"   --max_paths 5   --max_steps 7   --max_structures 150   --min_diameter 8   --max_diameter 8
```

---

## Implementation

- Builds a **graph of metabolites** using curated KEGG pathways via `build_graph_from_tsv()`.
- Standardises input InChIs to InChIKeys for unambiguous node representation (`Std.standardize`).
- Computes k-shortest paths via Dijkstra-like traversal (`k_shortest_paths`).
- If the start compound is absent from the network, calls **RetroPath2.0** through `RP.django_call()` to generate new reactions.
- Automatically merges RetroPath edges into the graph and retries path search.
- Produces a `.csv` summary of all retrieved paths with associated weights and compound identifiers.

### Key Arguments

- `--start_inchi`, `--end_inchi`: required; define the endpoints of the metabolic search. Our internal standardization pipeline protonates all molecules and collapses tautomers, in order to reduce the possibilities of not finding compounds in KEGG network and to allow the integration of RetroPath2.0 nodes. However, this might fail sometimes, and user is advised to try different InChI codes for the same molecule if nodes are not found in KEGG network.
  - example: pyruvate (InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)/p-1) might not be recognised as part of KEGG network, while pyruvic acid (InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)) might be.
- `--max_paths`: maximum number of shortest paths to return (default 3). Minimum 1.
- `--max_steps`: maximum RetroPath iteration depth (default 5). After 4 or 5 steps, memory and time usage increase exponentially. Minimum 1.
- `--max_structures`: maximum number of structures generated per step (default 100).
- `--min_diameter`, `--max_diameter`: control reaction rule diameter range (default 6–16). Minimum 6, maximum 16. The higher the diameter, the higher is the specificity of the reaction rule. Usually, a maximum diameter comprised between 6 and 8 generates more results but less specific for the substrate, while diameters above 12 generate less results but more specific.
- `--timeout`: execution limit in minutes for RetroPath subprocess.
- `--sink_file`, `--rules_file`: optional custom paths to RetroPath input data. See RetroPath2.0 original publication for more details on these. In `cape_app/algorithms/pathways/data`, a set of pre-parsed reaction rules and optimized sinks (`gsm_gapfilled_sink.csv`, `biocyc_cmpd_rpcleaned.csv`, `kegg_rp.csv`) for this software are provided.
  - `gsm_gapfilled_sink.csv` is a sink obtained by <i>Rhodococcus opacus</i> PD630 Genome Scale Model, with dead-end metabolites being removed. It is the most suitable sink for most applications.
  - `biocyc_cmpd_rpcleaned.csv` is a sink obtained by quering BioCyc for <i>R. opacus</i> metabolites (PD630 was not available). It is useful for increasing the possibility of merging RetroPath2.0 output with the KEGG network. Even though in most cases it would not be relevant, it is worth noting that not all compounds here are guaranteed to be actually present in PD630, and some compounds in PD630 might not be present in the sink.
  - `kegg_rp.csv` it is a sink obtained by considering all compounds in our KEGG network, including those of "white" edges and therefore not present in <i>R. opacus</i> PD630. When no other solution can be found with other types of sink, this one can offer a way to dramatically increase the chance to merge the networks, as reaching a sink compound here intrinsically means to merge the graphs.
- `--outdir`: output directory; created automatically if missing. <b>WARNING</b>: the content of this directory will be cleared before running the pathway prediction module.

---

## Example Output

pathway_id | pathway_weight | step_number | total_steps | source_compound | target_compound | EC_numbers| enzymes | enzyme_colors | reactions
|---------|-------------|------|---------|-------|---------|----------|----------------|---------|------|
1 | 16.213642968099766 | 1 | 7 | Terephthalate | "(3S,4R)-3,4-Dihydroxycyclohexa-1,5-diene-1,4-dicarboxylate" | 1.14.12.15; 1.18.1.- | "terephthalate 1,2-dioxygenase oxygenase component alpha subunit; terephthalate 1,2-dioxygenase reductase component"| white; white | "R05148, R05148, R05148, R05148, R05148, R05148; R05148, R05148, R05148"
1 | 16.213642968099766 | 2 | 7 | "(3S,4R)-3,4-Dihydroxycyclohexa-1,5-diene-1,4-dicarboxylate" | "3,4-Dihydroxybenzoate" | 1.3.1.53 | "1,2-dihydroxy-3,5-cyclohexadiene-1,4-dicarboxylate dehydrogenase" | white | "R01633, R01633, R01633"
1 | 16.213642968099766 | 3 | 7 | "3,4-Dihydroxybenzoate" | "3-Carboxy-cis,cis-muconate" | 1.13.11.3 | "protocatechuate 3,4-dioxygenase, alpha subunit" | green | "R01631, R01631, R01631, R01631, R01631, R01631, R01631, R01631"



Pathways aved to: `cape_app/pipeline_results/rp_results/pathways.csv`
Found 3 pathway(s) with 24 total steps

---

## Customisation Tips

- **Different rule sets**: swap `--rules_file` with any RetroRules export compatible with RetroPath2.0 to target specific reaction classes.
- **Update KEGG networks**: you can upload new tsv files inside `cape_app/algorithms/pathways/data/Pathways`, provided that they contain the following fields: enzyme_name, EC, source_inchi, source_inchi_key, target_inchi, target_inchi_key, source_name, target_name, reaction, color. The name of the file(s) must end with edges.tsv. Then run `./cape_app/algorithms/pathways/data/update_edge_weights.py`. Be careful: if run incorrectly, the folder might become corrupted, and consequently the software might not work properly anymore.
- **Heavier exploration**: increase `--max_structures` and `--max_steps` for complex transformations (at the cost of runtime), and/or lower `max_diameter`.
- **Batch usage**: combine multiple runs in shell loops or workflow managers.

---

[← Back to Home](../README.md)
