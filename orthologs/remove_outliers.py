import statistics
from fetchers.ncbi_filtering import fetch_protein_ids, fetch_sequences, parse_fasta

def infer_expected_length(sequences_dict):
    lengths = [len(seq) for seq in sequences_dict.values()]
    try:
        expected = statistics.mode(lengths)
    except statistics.StatisticsError:
        expected = statistics.median(lengths)
    return expected
'''
def length_score(actual_len, expected_len):
    if actual_len < 0.6 * expected_len:
        return -10  # likely fragment
    elif actual_len > 1.5 * expected_len:
        return -5   # likely fusion
    elif abs(actual_len - expected_len) < 0.15 * expected_len:
        return +5   # ideal
    else:
        return 0
'''
    
def is_suspect_header(header):
    """Check if a header contains keywords indicating low confidence annotations."""
    lower = header.lower()
    return any(keyword in lower for keyword in ['hypothetical', 'fragment', 'partial', 'putative'])   #can be extended

def iqr_bounds(lengths, k=1.5):
    """Return (lower, upper) bounds for outlier detection using Tukey's IQR rule."""
    q1, q3 = statistics.quantiles(lengths, n=4)[0], statistics.quantiles(lengths, n=4)[2]
    iqr = q3 - q1
    return q1 - k * iqr, q3 + k * iqr

def filter_by_length(sequences_dict):
    """Remove obvious length outliers using a robust IQR rule.
    * Records flagged as hypothetical/fragment/partial/putative are discarded first.
    * If <3 sequences remain, no statistical filtering is possible - return them.
    * If ≥5 sequences remain, use Tukey's rule (1.5 x IQR) to define bounds.
    * If 3-4 sequences remain, fall back to a broad window [0.5x, 2x] around the expected length (mode or median).
    """
    # --- step 1: drop low‑confidence annotations -----------------------------
    clean_dict = {
        header: seq
        for header, seq in sequences_dict.items()
        if not is_suspect_header(header)
    }

    if len(clean_dict) < 3:
        print("Too few sequences to apply length filter. Returning all non-suspect sequences.")
        return clean_dict

    # --- step 2: determine statistical bounds --------------------------------
    lengths = [len(seq) for seq in clean_dict.values()]

    if len(lengths) >= 5:
        lower, upper = iqr_bounds(lengths)
        print(f"Applying IQR filter: {lower:.1f} ≤ len ≤ {upper:.1f}")
    else:
        expected_len = infer_expected_length(clean_dict)
        lower, upper = 0.5 * expected_len, 2.0 * expected_len
        print(f"Small-sample fallback window: {lower:.1f} ≤ len ≤ {upper:.1f}")

    # --- step 3: filter -------------------------------------------------------
    filtered = {}
    for header, seq in clean_dict.items():
        l = len(seq)
        if lower <= l <= upper:
            filtered[header] = seq
        else:
            print(f"Removed length outlier: {header[:50]}... | Length: {l}")

    return filtered

if __name__ == '__main__':
    ec = '1.14.12.3'
    org = 'Actinomycetes'
    ec = '6.2.1.25'
    org = 'Nocardiaceae'
    protein_ids = fetch_protein_ids(ec, org)
    if not protein_ids:
        print("Nessun ID proteina trovato.")
        exit()
    fasta = fetch_sequences(protein_ids)
    parsed_sequences = parse_fasta(fasta)

    filtered_sequences = filter_by_length(parsed_sequences)
    print(f"\nRemaining sequences after outlier removal: {len(filtered_sequences)}")
    for header, seq in filtered_sequences.items():
        print(f">{header}\n{seq}\n")