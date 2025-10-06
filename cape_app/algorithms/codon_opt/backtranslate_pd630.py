"""
Back-translate protein FASTA to codon-optimized DNA for Rhodococcus opacus PD630.

Features
- Deterministic (max-frequency) or stochastic (frequency-weighted) codon choice
- Avoid internal restriction sites (e.g., TypeIIS RFC1000: BsaI/BsmBI/BbsI) on either strand
- Clean FASTA parser/writer; preserves headers
- Extensible for other hosts by swapping the CODON_FREQ table

Codon usage source: Kazusa codon usage database for 'species=37919' (R. opacus) - per-thousand frequencies.
See: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=37919
"""

import argparse
import random
import sys
from typing import Dict, List, Tuple, Iterable, Optional
from pathlib import Path
import os

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
    # Type IIS RFC1000
    # BsaI: GGTCTC, SapI: GCTCTTC
    "typeIIS": ["GGTCTC", "GCTCTTC"],
    # BioBrick RFC10
    # EcoRI: GAATTC, XbaI: TCTAGA, SpeI: ACTAGT, PstI: CTGCAG, NotI: GCGGCCGC
    "biobrick": ["GAATTC", "TCTAGA", "ACTAGT", "CTGCAG", "GCGGCCGC"],
    # No preset (user provides custom via --forbid)
    "none": [],
}

# Default restriction sites (kept for backward compatibility if user does not specify --assembly)
DEFAULT_FORBIDDEN = ["GGTCTC", "GCTCTTC"]  # BsaI, SapI (RFC1000 core)


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
        # Stochastic weighted choice: shuffle and sort by random weighted value
        codon_list, weights = zip(*codon_weights)

        # Sample all codons in random order weighted by frequency
        result = []
        available = list(codon_list)
        available_weights = list(weights)

        while available:
            total = sum(available_weights)
            probabilities = [w / total for w in available_weights]

            # Choose one codon based on weights
            chosen_idx = rng.choices(range(len(available)), weights=probabilities)[0]
            result.append(available[chosen_idx])

            # Remove chosen codon from available list
            available.pop(chosen_idx)
            available_weights.pop(chosen_idx)

        return result
    
    else:
        raise ValueError("Mode must be 'max' or 'weighted'")


def creates_forbidden_site(new_codon: str, dna_sequence: str, blocked_sites: List[str]) -> bool:
    """
    Check if sequence contains any forbidden restriction sites.
    
    Args:
        new_codon: newly added codon (string)
        dna_sequence: DNA sequence to check (string)
        blocked_sites: List of forbidden sites (strings)
    Returns:
        True if any forbidden site is found, False otherwise
    """
    if not blocked_sites:
        return False
    
    test_sequence = dna_sequence + new_codon

    return any(site in test_sequence for site in blocked_sites)


def backtranslate_protein(protein: str, mode: str, rng: random.Random = None, 
                         forbidden_sites: List[str] = None, max_backtrack: int = 10) -> str:
    """
    Convert protein sequence to DNA using codon optimization.
    Avoids creating forbidden restriction sites.

    Args:
        protein: Protein sequence (string; single-letter amino acids, '*' for stop)
        mode: "max" for deterministic max-frequency, "weighted" for stochastic weighted choice
        rng: Random number generator instance
        forbidden_sites: List of restriction sites to avoid (strings)
        max_backtrack: Maximum backtracking attempts per amino acid (default: 10)
    Returns:
        Back-translated DNA sequence (string)
    Raises:
        ValueError: If internal stop codons are found or invalid amino acids are present
        RuntimeError: If no valid codon can be placed
    """
    if rng is None:
        rng = random.Random(1)
    if forbidden_sites is None:
        forbidden_sites = DEFAULT_FORBIDDEN.copy()

    blocked_sites = get_blocked_sites(forbidden_sites)
    dna_parts = []
    
    # Clean protein sequence
    protein = protein.strip().upper()
    
    # Check for internal stop codons
    if "*" in protein[:-1]:
        raise ValueError("Internal stop codon '*' found in protein sequence")
    
    # Check presence of stop codon at end
    if protein[-1] != "*":
        protein += "*"
    
    # Track stats
    total_retries = 0
    total_backtracks = 0
    max_retries_for_position = 0

    # Stack-based backtracking: store (position, codon_choices, current_index)
    stack: List[List] = []
    position = 0
    while position < len(protein):
        aa = protein[position]

        # Get codon choices for this position (or from stack if backtracking)
        if position < len(stack):
            _, codon_choices, last_choice_idx = stack[position]
            start_idx = last_choice_idx + 1 # Try next codon
        else:
            # New position, get codon choices
            if aa == "*":
                codon_choices = choose_best_codons("*", mode, rng)
            else:
                codon_choices = choose_best_codons(aa, mode, rng)
            start_idx = 0
            stack.append([position, codon_choices, -1]) # -1 means none chosen yet

        # Try to place a codon
        codon_placed = False
        current_dna = "".join([stack[i][1][stack[i][2]] for i in range(position)])

        for choice_idx in range(start_idx, len(codon_choices)):
            codon = codon_choices[choice_idx]
            
            if not creates_forbidden_site(codon, current_dna, blocked_sites):
                # Success! Place this codon
                stack[position][2] = choice_idx  # Update choice index
                codon_placed = True
                
                if choice_idx > 0:
                    #print(f"    Position {position+1} ({aa}): used codon attempt {choice_idx+1}/{len(codon_choices)}: {codon}")
                    total_retries += choice_idx
                    max_retries_for_position = max(max_retries_for_position, choice_idx)
                
                position += 1
                break
        
        if not codon_placed:
            # No valid codon found - need to backtrack
            if position == 0:
                raise RuntimeError(
                    f"Cannot place any codon for first amino acid '{aa}'. "
                    f"All synonymous codons create forbidden sites. "
                    f"Try relaxing restriction site constraints."
                )
            
            # Calculate how far back we can go
            backtrack_distance = min(position, max_backtrack)
            backtrack_to = position - backtrack_distance
            
            #print(f"    ⚠ Position {position+1} ({aa}): all {len(codon_choices)} codons blocked, backtracking {backtrack_distance} positions to {backtrack_to+1}")
            total_backtracks += 1
            
            # Remove positions we're backtracking from
            stack = stack[:backtrack_to + 1]
            position = backtrack_to
            
            # Check if we've exhausted all options at the backtrack position
            if position >= 0 and stack[position][2] >= len(stack[position][1]) - 1:
                # All options exhausted at backtrack position too - need to backtrack further
                #print(f"    ⚠ Position {position+1} also exhausted, continuing backtrack...")
                continue
    
    # Build final sequence
    final_dna = "".join([stack[i][1][stack[i][2]] for i in range(len(protein))])
    
    #gc_content = (final_dna.count('G') + final_dna.count('C')) / len(final_dna) * 100
    
    # print(f"  ✓ Translation complete:")
    # print(f"    - DNA length: {len(final_dna)} bp ({len(final_dna)//3} codons)")
    # print(f"    - GC content: {gc_content:.1f}%")
    # print(f"    - Total codon retries: {total_retries}")
    # print(f"    - Total backtracks: {total_backtracks}")
    # print(f"    - Max retries at single position: {max_retries_for_position}")
    
    # Verify no forbidden sites in final sequence
    #found_sites = [site for site in blocked_sites if site in final_dna]
    #if found_sites:
        #print(f"  ⚠ WARNING: Found forbidden sites in final sequence: {found_sites}")
    #else:
        #print(f"  ✓ No forbidden restriction sites found")
    
    return final_dna


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
    parent_directory = output_file.parent
    os.makedirs(parent_directory, exist_ok=True)
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
                       help="Restriction sites to avoid (default: Type IIS RFC1000 sites)")
    parser.add_argument("--max-backtrack", type=int, default=10,
                   help="Maximum positions to backtrack when stuck (default: 10)")
    parser.add_argument(
        "--assembly",
        choices=list(ASSEMBLY_PRESETS.keys()),
        nargs="+",
        default=["none"],
        help=(
            "Assembly preset(s) to apply for illegal site exclusion. "
            "Options: typeIIS (BsaI,SapI), biobrick (EcoRI,XbaI,SpeI,PstI,NotI), none. "
            "You can pass multiple, e.g. --assembly typeIIS biobrick"
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
            dna_seq = backtranslate_protein(protein_seq, args.mode, rng, combined_sites, args.max_backtrack)
            dna_entries.append((header, dna_seq))
        
        # Write output
        write_fasta(dna_entries, args.output)
        print(f"DNA sequences written to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()