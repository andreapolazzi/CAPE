#!/usr/bin/env python3
"""
Run SelenzymeRF (selenzy.py) from Python and get a pandas DataFrame back.

Requirements (in your conda env):
  - pandas
  - the SelenzymeRF stack you already installed (rdkit, rxnmapper, etc.)

Example:
    df = run_selenzyme(
        smarts="[H]Oc1c([H])c(OC([H])([H])[H])c([H])c([H])c1C(=O)c1c([H])c([H])c([H])c([H])c1[H]>>[H]OC(=O)c1c([H])c([H])c(C(=O)c2c([H])c([H])c(OC([H])([H])[H])c([H])c2O[H])c([H])c1[H].[H]O[H]",
        repo_root="/work/igem/selenzyme2023",
        host_taxid=37919,                      # optional
        top_n=50,                              # optional
        outfile="selenzyme_results.csv"        # optional
    )
    print(df.head())
"""

from __future__ import annotations
import os, sys, subprocess
from pathlib import Path
import pandas as pd
from typing import Optional

def run_selenzyme(
    smarts: str,
    repo_root: str,
    datadir: Optional[str] = None,
    outdir: Optional[str] = None,
    host_taxid: Optional[int] = None,
    top_n: int = 50,
    nomsa: bool = True,
    outfile: str = "selenzyme_results.csv",
    python_exe: Optional[str] = None,
    timeout: Optional[int] = None,
) -> pd.DataFrame:
    """
    Execute selenzy.py with a SMARTS reaction and return the results as a DataFrame.
    """
    # Locate selenzy.py
    repo_root = Path(repo_root).resolve()
    selenzy = repo_root / "gitcode2023" / "selenzyPro" / "selenzy.py"
    if not selenzy.is_file():
        raise FileNotFoundError(f"Could not find selenzy.py at: {selenzy}")

    # Set up data and output directories
    if datadir is None:
        datadir = repo_root / "gitcode2023" / "selenzyPro" / "data" / "data_2023"
    else:
        datadir = Path(datadir)
        
    if outdir is None:
        outdir = repo_root / "runout"
    else:
        outdir = Path(outdir)

    # Create outdir if needed
    outdir.mkdir(parents=True, exist_ok=True)    

    # Ensure outdir has a trailing separator (some tooling expects it)
    outdir_with_sep = str(outdir) + os.sep

    # Build command
    py = python_exe or sys.executable
    cmd = [py, str(selenzy),
           "-smarts", smarts,
           str(datadir),
           outdir_with_sep,
           "-tar", str(top_n)]
    if nomsa:
        cmd += ["-NoMSA"]
    if host_taxid is not None:
        cmd += ["-host", str(host_taxid)]
    cmd += ["-outfile", outfile]

    # Run it
    proc = subprocess.run(
        cmd, cwd=selenzy.parent, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=timeout
    )

    # Figure out which CSV to read
    csv_path = outdir / outfile
    if not csv_path.exists():
        # Fallback: pick the newest CSV in outdir
        candidates = sorted(outdir.glob("*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
        if candidates:
            csv_path = candidates[0]
        else:
            raise FileNotFoundError(f"Expected CSV not found: {outdir/outfile}")

    # Load DataFrame
    df = pd.read_csv(csv_path)

    return df



if __name__ == "__main__":
    # Minimal CLI usage:
    import argparse
    p = argparse.ArgumentParser(description="Run SelenzymeRF and print top rows.")
    p.add_argument("--smarts", required=True, help="Reaction SMARTS/SMIRKS")
    p.add_argument("--repo-root", default=".")
    p.add_argument("--datadir", default=None)
    p.add_argument("--outdir", default=None)
    p.add_argument("--host", type=int, default=None)
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--outfile", default="selenzyme_results.csv")

    # MSA is disabled by default; allow explicit enable/disable
    g = p.add_mutually_exclusive_group()
    g.add_argument("--msa", dest="nomsa", action="store_false", help="Enable MSA (default: disabled)")
    g.add_argument("--no-msa", dest="nomsa", action="store_true", help="Disable MSA (default)")
    p.set_defaults(nomsa=True)
    args = p.parse_args()

    df = run_selenzyme(
        smarts=args.smarts,
        repo_root=args.repo_root,
        datadir=args.datadir,
        outdir=args.outdir,
        host_taxid=args.host,
        top_n=args.top,
        nomsa=args.nomsa,
        outfile=args.outfile
    )
    # Report the expected CSV path
    from pathlib import Path as _P
    _repo = _P(args.repo_root).resolve()
    _out = _P(args.outdir) if args.outdir else (_repo / "runout")
    _csv = _out / args.outfile
    print(f"[OK] CSV: {_csv}")
    # Show a quick preview
    show_cols = [c for c in df.columns if any(k in c.lower() for k in ["ec","uniprot","organism","score","rxn sim","tax"])]
    print(df[show_cols or df.columns].head(20).to_string(index=False))
