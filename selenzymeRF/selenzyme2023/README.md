# SelenzymeRF Runner

A Python wrapper (`run_selenzyme.py`) to run **SelenzymeRF** and obtain results in `pandas.DataFrame` or CSV format

---

## 1. Setup environment

You must have [conda] installed.
Then, execute the setup script:

```bash
bash setup_env.sh
```

conda activate selenzyme

---

## 2. Run from CLI

You can run the tool via the provided script:

```bash
python run_selenzyme.py \
  --smarts "C=C>>CC" \
  --repo-root . \
  --host 37919 \
  --top 20 \
  --outfile results.csv
```

### Arguments

- `--smarts` (**required**) : SMARTS/SMIRKS of the reaction
- `--repo-root` : path to the repo root (default: `.`)
- `--datadir` : data directory (default: `<repo>/gitcode2023/selenzyPro/data/data_2023`)
- `--outdir` : output directory (default: `<repo>/runout`)
- `--host` : NCBI TaxID for host organism (default: `37919`)
- `--top` : number of top results to keep (default: `50`)
- `--outfile` : name of the CSV output file (default: `selenzyme_results.csv`)
- `--msa` / `--no-msa` : enable or disable MSA (disabled by default)

---

## 3. Output

The script will:
- Print the top rows of the results to the console
- Save a CSV file (`--outfile`) in the specified `--outdir`

Example:

```
[OK] CSV written: ./runout/results.csv
   ec_number   uniprot   organism   score
   1.1.1.1     P12345    E.coli     0.87
   ...
```

---

## 4. Use as Python library

You can also call the function directly in your own code:

```python
from run_selenzyme import run_selenzyme

df = run_selenzyme(
    smarts="C=C>>CC",
    repo_root=".",
    host_taxid=37919,
    top_n=20
)

print(df.head())
```