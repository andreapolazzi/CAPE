#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Back-translate protein FASTA to codon-optimized DNA for Rhodococcus opacus PD630.

Features
- Deterministic (max-frequency) or stochastic (frequency-weighted) codon choice
- Avoid internal restriction sites (e.g., Golden Gate: BsaI/BsmBI/BbsI) on either strand
- Clean FASTA parser/writer; preserves headers
- Extensible for other hosts by swapping the CODON_FREQ table

Codon usage source: Kazusa codon usage database for 'species=37919' (R. opacus) – per-thousand frequencies.
See: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37919
"""

import argparse
import random
import sys
from typing import Dict, List, Tuple, Iterable

# -------------------------
# Kazusa per-thousand frequencies for R. opacus (species=37919), RNA->DNA (U->T).
# Rows copied from the table; frequencies rounded as shown by Kazusa.
# -------------------------
CODON_FREQ_PD630: Dict[str, float] = {
    # U-start (T-start)
    "TTT": 2.9,  "TCT": 2.6,  "TAT": 4.3,  "TGT": 2.4,
    "TTC": 29.8, "TCC": 13.6, "TAC": 19.6, "TGC": 7.0,
    "TTA": 1.0,  "TCA": 4.0,  "TAA": 0.4,  "TGA": 2.6,
    "TTG": 8.3,  "TCG": 18.6, "TAG": 0.6,  "TGG": 15.2,
    # C-start
    "CTT": 5.7,  "CCT": 3.8,  "CAT": 4.9,  "CGT": 7.8,
    "CTC": 37.4, "CCC": 18.1, "CAC": 20.1, "CGC": 31.2,
    "CTA": 2.5,  "CCA": 5.4,  "CAA": 5.6,  "CGA": 8.3,
    "CTG": 38.5, "CCG": 29.3, "CAG": 26.0, "CGG": 18.4,
    # A-start
    "ATT": 6.1,  "ACT": 4.7,  "AAT": 5.5,  "AGT": 4.6,
    "ATC": 35.6, "ACC": 33.7, "AAC": 20.7, "AGC": 12.9,
    "ATA": 2.3,  "ACA": 6.9,  "AAA": 5.4,  "AGA": 1.7,
    "ATG": 21.5, "ACG": 17.8, "AAG": 21.8, "AGG": 3.4,
    # G-start
    "GTT": 7.9,  "GCT": 9.2,  "GAT": 15.3, "GGT": 17.4,
    "GTC": 38.3, "GCC": 45.9, "GAC": 49.0, "GGC": 41.3,
    "GTA": 3.6,  "GCA": 15.0, "GAA": 20.2, "GGA": 12.6,
    "GTG": 29.1, "GCG": 37.4, "GAG": 41.8, "GGG": 15.5,
}

# Standard genetic code (DNA codons -> amino acid one-letter; * for stop).
# (NCBI code 11 – the usual “bacterial”/standard code for nuclear genes.)
CODON_TO_AA: Dict[str, str] = {
    # F/L
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # I/M/V
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # S/P/T/A
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Y/*/C/W
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGG":"W",
    # H/Q/N/K/D/E
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    # R/G
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def invert_code(codon_to_aa: Dict[str, str]) -> Dict[str, List[str]]:
    aa_to_codons: Dict[str, List[str]] = {}
    for codon, aa in codon_to_aa.items():
        aa_to_codons.setdefault(aa, []).append(codon)
    return aa_to_codons

AA_TO_CODONS = invert_code(CODON_TO_AA)

# Default forbidden sites (Golden Gate workhorses). Add more via CLI.
DEFAULT_FORBIDDEN = ["GGTCTC", "CGTCTC", "GAAGAC"]  # BsaI, BsmBI, BbsI recognition sites (sense)
# We'll automatically also block reverse complements.

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

def all_blocked_sites(sites: Iterable[str]) -> List[str]:
    s = []
    for site in sites:
        site = site.upper()
        s.append(site)
        s.append(revcomp(site))
    # Deduplicate while preserving order
    out, seen = [], set()
    for x in s:
        if x not in seen:
            out.append(x); seen.add(x)
    return out

def parse_fasta(fp) -> Iterable[Tuple[str, str]]:
    header, seq_chunks = None, []
    for line in fp:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.replace(" ", "").upper())
    if header is not None:
        yield header, "".join(seq_chunks)

def wrap(seq: str, width: int = 70) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def choose_codon_list(aa: str, mode: str, rng: random.Random) -> List[str]:
    """Return codons for aa sorted by preference (or a randomized weighted single choice)."""
    codons = AA_TO_CODONS.get(aa)
    if not codons:
        raise ValueError(f"Unsupported amino acid symbol: {aa!r}")
    # Use embedded frequencies (fallback to small epsilon if missing)
    items = [(c, CODON_FREQ_PD630.get(c, 0.0001)) for c in codons]
    if mode == "max":
        # Highest-frequency first
        return [c for c, _ in sorted(items, key=lambda x: x[1], reverse=True)]
    elif mode == "weighted":
        # Weighted random choice -> but we still return a preference list where the first is sampled
        codon_names, weights = zip(*items)
        total = sum(weights)
        probs = [w/total if total > 0 else 1.0/len(weights) for w in weights]
        # sample a top codon; then put remaining in descending order to use as fallbacks
        chosen = rng.choices(codon_names, weights=probs, k=1)[0]
        rest = [c for c, _ in sorted(items, key=lambda x: x[1], reverse=True) if c != chosen]
        return [chosen] + rest
    else:
        raise ValueError("mode must be 'max' or 'weighted'")

def would_create_forbidden(new_seq_tail: str, blocked_sites: List[str]) -> bool:
    s = new_seq_tail.upper()
    for site in blocked_sites:
        if site in s:
            return True
    return False

def backtranslate(protein: str, mode: str, rng: random.Random, forbidden_sites: List[str]) -> str:
    """
    Greedy builder: at each AA, try codons in preference order; skip any that
    would introduce a forbidden site in the *current* sequence (scan last 100 nt window).
    """
    blocked = all_blocked_sites(forbidden_sites)
    dna: List[str] = []
    # If trailing '*' present, keep as explicit STOP at the end (choose most common stop)
    stop_pref = choose_codon_list("*", mode="max", rng=rng)  # use deterministic for stop
    # Expect only one AA '*' mapping
    # (Kazusa stop usage in PD630: TGA more frequent than TAA/TAG)

    prot = protein.strip().upper()
    if "*" in prot[:-1]:
        raise ValueError("Internal stop '*' found before the end of the protein sequence.")

    for i, aa in enumerate(prot):
        if aa == "*":
            # End: append a STOP using the most frequent among TGA/TAA/TAG (already ordered)
            for codon in stop_pref:
                candidate = "".join(dna) + codon
                tail = candidate[-100:]
                if not would_create_forbidden(tail, blocked):
                    dna.append(codon)
                    break
            else:
                # If all stops conflict, just use the top stop anyway
                dna.append(stop_pref[0])
            continue

        # Try preferred codons for this AA
        tried = choose_codon_list(aa, mode, rng)
        placed = False
        for codon in tried:
            candidate = "".join(dna) + codon
            tail = candidate[-100:]
            if not would_create_forbidden(tail, blocked):
                dna.append(codon)
                placed = True
                break
        if not placed:
            # Last-resort: place the top codon even if it makes a site; user can resolve later
            dna.append(tried[0])

    return "".join(dna)

def main():
    ap = argparse.ArgumentParser(
        description="Back-translate protein FASTA to PD630-codon-optimized DNA with optional restriction-site exclusion."
    )
    ap.add_argument("-i", "--input", required=True, help="Protein FASTA file")
    ap.add_argument("-o", "--output", required=True, help="Output DNA FASTA file")
    ap.add_argument("--mode", choices=["max", "weighted"], default="max",
                    help="Codon selection: 'max' (most frequent) or 'weighted' (sampled by frequency). Default: max.")
    ap.add_argument("--seed", type=int, default=1, help="Random seed (for weighted mode).")
    ap.add_argument("--forbid", nargs="*", default=DEFAULT_FORBIDDEN,
                    help="Recognition sequences to forbid internally (e.g. GGTCTC CGTCTC GAAGAC). "
                         "Reverse complements are automatically forbidden too.")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    with open(args.input) as fp:
        entries = list(parse_fasta(fp))

    out_lines = []
    for header, prot in entries:
        # Strip terminal stop if absent; allow adding one by including '*' in the protein
        dna = backtranslate(prot, mode=args.mode, rng=rng, forbidden_sites=args.forbid)
        out_lines.append(f">{header}")
        out_lines.append(wrap(dna))

    with open(args.output, "w") as out:
        out.write("\n".join(out_lines) + "\n")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)