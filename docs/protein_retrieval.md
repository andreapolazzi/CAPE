# Protein Retrieval Module

The orchestration script lives at:
`cape_app/algorithms/orthologs/processing.py`

Supporting NCBI wrappers:
`cape_app/algorithms/orthologs/ncbiprot_api.py`

This module queries the NCBI Protein E-utilities for enzyme sequences mapped to specific EC numbers, prioritises reviewed annotations in the Actinomycetota lineage, filters out low-confidence headers, and exports curated sequence tables for downstream screening.

The CLI depends on `requests`, `pandas`, and `regex` in addition to the standard library.

## Input

- One or more EC numbers supplied on the command line (e.g. `1.14.12.15`).
- Alternatively, a text file with one EC number per line (`--file path/to/ec_list.txt`). Lines starting with `#` are ignored.

## Output

- Console summary reporting the taxonomy level that first returned matches and preview tables of curated sequences (default behaviour).
- Optional CSVs written to `--outdir`, including an extra file for suspect entries when `--include-hidden` is set.

## Usage

```bash
# Process a single EC number, print summary + table preview
python cape_app/algorithms/orthologs/processing.py 1.14.12.15

# Process multiple EC numbers, retain only curated entries, and export to CSVs
python cape_app/algorithms/orthologs/processing.py \
  1.14.12.15 3.1.1.1 \
  --outdir results/protein_retrieval

# Read EC numbers from file and include suspect sequences in a separate CSV
python cape_app/algorithms/orthologs/processing.py \
  --file configs/ec_numbers.txt \
  --outdir results/protein_retrieval \
  --include-hidden
```

## Implementation

- Iteratively queries the NCBI Protein database for each EC number, starting from _Rhodococcus_ and widening through higher _Actinomycetota_ taxa until hits are found.
- Fetches sequence FASTA records in chunks, parses them into header/sequence pairs, and normalises whitespace.
- Flags suspect entries containing terms such as _hypothetical_, _fragment_, or _putative_; these entries are discarded unless `--include-hidden` is requested.
- Builds pandas DataFrames with organism names, cleaned headers, sequence lengths, and direct NCBI links for quick inspection.
- Writes CSV exports (main and optional hidden) and prints  progress updates to stdout.

### Key Arguments

- `ec_numbers`: positional list of EC numbers to retrieve when `--file` is not supplied.
- `--file` / `-f`: point to a file of EC numbers (one per line) to avoid long CLI invocations.
- `--outdir` / `-o`: directory in which to save CSV results; if omitted, results print to stdout only.
- `--include-hidden`: instructs the exporter to write suspect entries alongside the curated set.
- `batch_size` / `chunk_size`: tunable fetch sizes exposed via `ncbiprot_api.fetch_protein_ids` and `fetch_sequences` should you call the functions directly.

### Customisation Tips

- **Taxonomy breadth**: adjust the `orgs` list in `process_ec_number` to target different host clades or prioritise other organisms.
- **Header triage**: expand `is_suspect_header` with additional keywords relevant to your curation rules.
- **Rate limiting**: insert sleeps or API keys if you hit NCBI usage ceilings; `requests` session handling can also be added for retries.

[← Back to Home](../README.md)
