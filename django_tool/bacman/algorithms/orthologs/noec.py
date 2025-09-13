"""noec.py — minimal wrapper to call SelenzymeRF (selenzy.py) via run_selenzyme.py
and return a pandas DataFrame, with optional column selection and batch helpers.

Designed for use inside the Django backend. Keep it stateless and side-effect free.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Optional, Dict, Any, Union
import importlib.util
import pandas as pd


# --- Dynamic loader for run_selenzyme.py ---
_DEF_REPO_REL = Path(__file__).resolve().parents[4] / "selenzymeRF" / "selenzyme2023"
_DEF_RUN_PATH = _DEF_REPO_REL / "run_selenzyme.py"


def _load_run_selenzyme(run_path: Union[str, Path] = _DEF_RUN_PATH):
    """Load run_selenzyme.py from a file path and return the module.
    Allows this module to be used from different app working dirs.
    """
    run_path = Path(run_path)
    if not run_path.is_file():
        raise FileNotFoundError(f"run_selenzyme.py not found at: {run_path}")
    spec = importlib.util.spec_from_file_location("run_selenzyme", str(run_path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load spec for: {run_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


# --- Core single-query function ---

def run_noec(
    smarts: str,
    repo_root: Optional[Union[str, Path]] = None,
    *,
    datadir: Optional[Union[str, Path]] = None,
    outdir: Optional[Union[str, Path]] = None,
    host_taxid: Optional[int] = None,
    top_n: int = 50,
    nomsa: bool = True,
    outfile: str = "selenzyme_results.csv",
    python_exe: Optional[str] = None,
    timeout: Optional[int] = None,
    include: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Run SelenzymeRF for a single SMARTS string and return a DataFrame.

    If *include* is given, select only those columns (silently ignoring missing).

    All other kwargs are forwarded to run_selenzyme.run_selenzyme.
    """
    mod = _load_run_selenzyme()
    repo_root = Path(repo_root) if repo_root is not None else str(_DEF_REPO_REL)

    res = mod.run_selenzyme(
        smarts=smarts,
        repo_root=str(repo_root),
        datadir=str(datadir) if datadir else None,
        outdir=str(outdir) if outdir else None,
        host_taxid=host_taxid,
        top_n=top_n,
        nomsa=nomsa,
        outfile=outfile,
        python_exe=python_exe,
        timeout=timeout,
    )

    # Expect a DataFrame from run_selenzyme
    df = res
    assert isinstance(df, pd.DataFrame)

    if include is not None:
        wanted = [c for c in include if c in df.columns]
        if wanted:
            df = df[wanted]

    return df


# --- Batch helpers ---

def run_noec_batch(
    smarts_list: Iterable[str],
    repo_root: Optional[Union[str, Path]] = None,
    *,
    include: Optional[Iterable[str]] = None,
    **kwargs: Any,
) -> Dict[str, pd.DataFrame]:
    """Process multiple SMARTS sequentially.

    Returns a dict mapping each SMARTS to its DataFrame.
    """
    results: Dict[str, pd.DataFrame] = {}
    for s in smarts_list:
        df = run_noec(
            s,
            repo_root=repo_root,
            include=include,
            **kwargs,
        )
        results[s] = df

    return results

if __name__ == "__main__":
    # All columns : ["Seq. ID","Description","Organism Source","Tax. distance","Rxn. ID","Rxn. ID - additional","EC Number","Uniprot protein evidence","Consv. Score","Rxn Sim.","Rxn Sim RF.","Direction Used","Combined score","% helices","% sheets","% turns","% coils","Mol. Weight","Isoelec. Point","Polar %","Query","Hit"]

    # test_smarts = "[H]Oc1c([H])c(OC([H])([H])[H])c([H])c([H])c1C(=O)c1c([H])c([H])c([H])c([H])c1[H]>>[H]OC(=O)c1c([H])c([H])c(C(=O)c2c([H])c([H])c(OC([H])([H])[H])c([H])c2O[H])c([H])c1[H].[H]O[H]"
    # df = run_noec(
    #     smarts=test_smarts,
    #     repo_root=_DEF_REPO_REL,
    #     host_taxid=37919,
    #     include=["Description", "Organism Source", "Tax. distance", "EC Number", "Consv. Score", "Rxn Sim.", "Combined score"],
    # )

    # print(df.head())

    test_multi = ["[H]Oc1c([H])c(OC([H])([H])[H])c([H])c([H])c1C(=O)c1c([H])c([H])c([H])c([H])c1[H]>>[H]OC(=O)c1c([H])c([H])c(C(=O)c2c([H])c([H])c(OC([H])([H])[H])c([H])c2O[H])c([H])c1[H].[H]O[H]",
                    "[H]C1=CC(=CC=C1C=O)O[H]>>[H]C1=CC(=CC=C1C(=O)O[H])O[H]"]
    
    res = run_noec_batch(
        test_multi,
        repo_root=_DEF_REPO_REL,
        host_taxid=37919,
        include=["Description", "Organism Source", "Tax. distance", "EC Number", "Consv. Score", "Rxn Sim.", "Combined score"],
    )

    for smarts, df in res.items():
        print(f"Results for SMARTS: {smarts}")
        print(df.head())