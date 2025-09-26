"""
Back-translate protein FASTA to codon-optimized DNA for Rhodococcus opacus PD630.

Features
- Deterministic (max-frequency) or stochastic (frequency-weighted) codon choice
- Avoid internal restriction sites (e.g., Golden Gate: BsaI/BsmBI/BbsI) on either strand
- Clean FASTA parser/writer; preserves headers
- Extensible for other hosts by swapping the CODON_FREQ table

Codon usage source: Kazusa codon usage database for 'species=37919' (R. opacus) - per-thousand frequencies.
See: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37919
"""

import argparse
import random
import sys
from typing import Dict, List, Tuple, Iterable, Optional

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

# Create reverse mapping (amino acid -> codons)
AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# Assembly presets: recognition sequences to forbid (sense; reverse complements will be auto-blocked)
ASSEMBLY_PRESETS: Dict[str, List[str]] = {
    # Type IIS Golden Gate
    # BsaI: GGTCTC, SapI: GCTCTTC
    "golden_gate": ["GGTCTC", "GCTCTTC"],
    # BioBrick RFC10
    # EcoRI: GAATTC, XbaI: TCTAGA, SpeI: ACTAGT, PstI: CTGCAG, NotI: GCGGCCGC
    "biobrick": ["GAATTC", "TCTAGA", "ACTAGT", "CTGCAG", "GCGGCCGC"],
    # No preset (user provides custom via --forbid)
    "none": [],
}

# Default restriction sites (kept for backward compatibility if user does not specify --assembly)
DEFAULT_FORBIDDEN = ["GGTCTC", "GCTCTTC"]  # BsaI, SapI (Golden Gate core)


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.
    
    Args:
        seq: Input DNA sequence (string)
    Returns:
        Reverse complement sequence (string)
    """
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(complement[base] for base in reversed(seq.upper()))


def get_blocked_sites(forbidden_sites: List[str]) -> List[str]:
    """
    Get list of sites to avoid including reverse complements; deduplicate preserving order.
    
    Args:
        forbidden_sites: List of restriction sites (strings)
    Returns:
        List of unique sites and their reverse complements (strings)
    """
    seen = set()
    blocked: List[str] = []
    for site in forbidden_sites:
        s = site.upper()
        rc = reverse_complement(s)
        for x in (s, rc):
            if x not in seen:
                blocked.append(x)
                seen.add(x)
    return blocked


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
    """
    Parse FASTA file and return list of (header, sequence) tuples.
    
    Args:
        filepath: Path to FASTA file
    Returns:
        List of (header, sequence) tuples
    """
    entries = []
    with open(filepath) as f:
        header, sequence = None, []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.replace(" ", "").upper())
        
        if header:
            entries.append((header, "".join(sequence)))
    
    return entries


def format_fasta_sequence(seq: str, width: int = 70) -> str:
    """
    Format sequence with line breaks.
    
    Args:
        seq: Input sequence (string)
        width: Line width (default: 70)
    Returns:
        Formatted sequence with line breaks (string)
    """
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def choose_best_codons(aa: str, mode: str, rng: random.Random) -> List[str]:
    """
    Return codons for amino acid ordered by preference.
    
    Args:
        aa: Amino acid (single-letter code; '*' for stop)
        mode: "max" for deterministic max-frequency, "weighted" for stochastic weighted choice
        rng: Random number generator instance
    Returns:
        List of codons ordered by preference
    Raises:
        ValueError: If amino acid is unknown or mode is invalid
    """
    if aa not in AA_TO_CODONS:
        raise ValueError(f"Unknown amino acid: {aa}")
    
    codons = AA_TO_CODONS[aa]
    codon_weights = [(codon, CODON_FREQ_PD630.get(codon, 0.0001)) for codon in codons]
    
    if mode == "max":
        # Sort by frequency (highest first)
        return [codon for codon, _ in sorted(codon_weights, key=lambda x: x[1], reverse=True)]
    
    elif mode == "weighted":
        # Random weighted choice, then sort rest by frequency
        codons_list, weights = zip(*codon_weights)
        total_weight = sum(weights)
        probabilities = [w / total_weight for w in weights]
        
        chosen = rng.choices(codons_list, weights=probabilities)[0]
        remaining = [codon for codon, _ in sorted(codon_weights, key=lambda x: x[1], reverse=True) 
                    if codon != chosen]
        
        return [chosen] + remaining
    
    else:
        raise ValueError("Mode must be 'max' or 'weighted'")


def creates_forbidden_site(dna_sequence: str, blocked_sites: List[str]) -> bool:
    """
    Check if sequence contains any forbidden restriction sites.
    
    Args:
        dna_sequence: DNA sequence to check (string)
        blocked_sites: List of forbidden sites (strings)
    Returns:
        True if any forbidden site is found, False otherwise
    Raises:
        ValueError: If dna_sequence contains invalid characters
    """
    seq = dna_sequence.upper()
    if any(base not in "ACGT" for base in seq):
        raise ValueError("DNA sequence contains invalid characters (only A,C,G,T allowed)")
    return any(site in seq for site in blocked_sites)


def backtranslate_protein(protein: str, mode: str, rng: random.Random, 
                         forbidden_sites: List[str]) -> str:
    """
    Convert protein sequence to DNA using codon optimization.
    Avoids creating forbidden restriction sites.

    Args:
        protein: Protein sequence (string; single-letter amino acids, '*' for stop)
        mode: "max" for deterministic max-frequency, "weighted" for stochastic weighted choice
        rng: Random number generator instance
        forbidden_sites: List of restriction sites to avoid (strings)
    Returns:
        Back-translated DNA sequence (string)
    Raises:
        ValueError: If internal stop codons are found or invalid amino acids are present
        RuntimeError: If no valid codon can be placed
    """

    blocked_sites = get_blocked_sites(forbidden_sites)
    dna_parts = []
    
    # Clean protein sequence
    protein = protein.strip().upper()
    
    # Check for internal stop codons
    if "*" in protein[:-1]:
        raise ValueError("Internal stop codon '*' found in protein sequence")
    
    for i, aa in enumerate(protein):
        if aa == "*":
            # Use most frequent stop codon
            preferred_codons = choose_best_codons("*", "max", rng)
        else:
            preferred_codons = choose_best_codons(aa, mode, rng)
        
        # Try each codon until we find one that doesn't create forbidden sites
        codon_placed = False
        for codon in preferred_codons:
            test_sequence = "".join(dna_parts) + codon
            # Check last 100 nucleotides for restriction sites
            check_region = test_sequence[-100:]
            
            if not creates_forbidden_site(check_region, blocked_sites):
                dna_parts.append(codon)
                codon_placed = True
                break
        
        if not codon_placed:
            raise RuntimeError(f"No valid codon could be placed for amino acid '{aa}' at position {i+1}. All synonyms create forbidden sites.")
    
    return "".join(dna_parts)


def write_fasta(entries: List[Tuple[str, str]], output_file: str):
    """
    Write sequences to FASTA file.
    
    Args:
        entries: List of (header, sequence) tuples
        output_file: Path to output FASTA file
    Returns:
        None
    Raises:
        IOError: If file cannot be written
    """
    with open(output_file, "w") as f:
        for header, sequence in entries:
            f.write(f">{header}\n")
            f.write(format_fasta_sequence(sequence) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Back-translate protein FASTA to PD630-codon-optimized DNA"
    )
    parser.add_argument("-i", "--input", required=True, 
                       help="Input protein FASTA file")
    parser.add_argument("-o", "--output", required=True, 
                       help="Output DNA FASTA file")
    parser.add_argument("--mode", choices=["max", "weighted"], default="max",
                       help="Codon selection mode (default: max)")
    parser.add_argument("--seed", type=int, default=1, 
                       help="Random seed for weighted mode (default: 1)")
    parser.add_argument("--forbid", nargs="*", default=DEFAULT_FORBIDDEN,
                       help="Restriction sites to avoid (default: Golden Gate sites)")
    parser.add_argument(
        "--assembly",
        choices=list(ASSEMBLY_PRESETS.keys()),
        nargs="+",
        default=["none"],
        help=(
            "Assembly preset(s) to apply for illegal site exclusion. "
            "Options: golden_gate (BsaI,SapI), biobrick (EcoRI,XbaI,SpeI,PstI,NotI), none. "
            "You can pass multiple, e.g. --assembly golden_gate biobrick"
        ),
    )
    
    args = parser.parse_args()
    
    try:
        rng = random.Random(args.seed)

        # Compute forbidden sites from assembly preset(s) and user-provided additions
        preset_sites: List[str] = []
        for preset in (args.assembly or []):
            preset_sites.extend(ASSEMBLY_PRESETS.get(preset, []))
        # Backward-compatibility: if user provided --forbid, merge them
        extra_sites = args.forbid or []
        combined_sites = preset_sites + list(extra_sites)
        
        # Parse input
        protein_entries = parse_fasta(args.input)
        print(f"Processing {len(protein_entries)} protein sequences...")
        
        # Back-translate each protein
        dna_entries = []
        for header, protein_seq in protein_entries:
            print(f"Back-translating: {header}")
            dna_seq = backtranslate_protein(protein_seq, args.mode, rng, combined_sites)
            dna_entries.append((header, dna_seq))
        
        # Write output
        write_fasta(dna_entries, args.output)
        print(f"DNA sequences written to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()