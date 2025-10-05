# Enzyme Selection Module

SelenzymeRF client classes and CLI live at:
`cape_app/algorithms/selenzyme/api_selenzyme.py`
`cape_app/algorithms/selenzyme/noec_api.py`

This module submits SMARTS transformations to the SelenzymeRF service, extracts enzyme recommendations ranked by sequence similarity and phylogenetic distance, and saves curated tables for candidate selection. It is network-bound and depends on `requests`, `pandas`, and the standard library.

## Input

- One or more reaction SMARTS strings passed directly on the command line.
- Alternatively, a text file containing SMARTS expressions (one per line, `#` for comments) supplied via `--file` or positional argument pointing to an existing file.

## Output

- Console summaries listing the top-scoring enzyme entries per SMARTS query (default).
- Optional CSV exports: single-file via `--output` for solitary SMARTS, or per-query multiple files written to `--outdir`.

## Usage

```bash
# Single SMARTS query; print top entries to the console
python cape_app/algorithms/selenzyme/noec_api.py "CC(=O)O>>CCO"

# Multiple SMARTS with CSV outputs (one file per query)
python cape_app/algorithms/selenzyme/noec_api.py \
  "CC(=O)O>>CCO" "c1ccccc1O>>c1ccccc1" \
  --targets 100 \
  --host 37919 \
  --outdir results/selenzymerf

# Load SMARTS from file and save the single result to a CSV
python cape_app/algorithms/selenzyme/noec_api.py \
  --file configs/smarts_list.txt \
  --output results/acylation_candidates.csv
```

## Implementation

- `SelenzymeRFClient` wraps the REST API: submits SMARTS payloads, handles timeouts, and parses the JSON field that bundles tabular data as nested dictionaries.
- The CLI normalises SMARTS input from files or positional arguments, issues rate-limited API calls (`--delay`) to respect service usage, and converts responses into pandas DataFrames.
- Helper utilities create per-query CSV exports, skip empty result sets, and surface warnings when requests fail.
- All console logging uses straightforward `print` statements so behaviour is transparent in pipelines.

### Key Arguments

- `smarts`: positional SMARTS patterns when no file is provided; recognised as a file path if the single argument exists on disk.
- `--file` / `-f`: path to a newline-delimited SMARTS list.
- `--targets` / `-t`: maximum number of enzyme hits to request (default 50).
- `--host`: NCBI taxonomic identifier for host weighting (default `37919`, *R. opacus*).
- `--timeout`: HTTP timeout in seconds; complements `requests` error handling.
- `--delay`: pause between API calls during batch runs to avoid rate-limit issues.
- `--output` / `--outdir`: control CSV export layout for single vs multiple SMARTS queries.
- `--api-url`: point to a different url (default <http://selenzymerf.synbiochem.co.uk>)

### Customisation Tips

- **Base URL overrides**: point `--api-url` to an alternate SelenzymeRF deployment or a local cache.
- **Retry strategy**: wrap `query_smarts` with backoff logic or swap to a persistent `requests.Session` for large-scale batch jobs.

[← Back to Home](../README.md)
