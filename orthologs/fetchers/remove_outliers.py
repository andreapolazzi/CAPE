import statistics
import numpy as np

def infer_expected_length(sequences_dict):
    """Infer the expected length of sequences based on the most common length.
    If there's no clear mode, use the median length instead.
    """
    lengths = [len(seq) for seq in sequences_dict.values()]
    try:
        expected = statistics.mode(lengths)
    except statistics.StatisticsError:
        expected = statistics.median(lengths)
    return expected
    
def is_suspect_header(header):
    """Check if a header contains keywords indicating low confidence annotations."""
    lower = header.lower()
    return any(keyword in lower for keyword in ['hypothetical', 'fragment', 'partial', 'putative'])   #can be extended

def iqr_bounds(lengths, k=1.5):
    """Return (lower, upper) bounds for outlier detection using Tukey's IQR rule, with fallback if IQR is 0 or bounds are illogical."""
    q1 = np.percentile(lengths, 25)
    q3 = np.percentile(lengths, 75)
    iqr = q3 - q1

    lower, upper = q1 - k * iqr, q3 + k * iqr

    # Fallback if IQR is 0, negative, or bounds are nonsensical
    if iqr <= 0 or lower < 0 or upper < 1.1 * q3:
        median = np.median(lengths)
        print("IQR unreliable or bounds too tight — using fallback window around median.")
        return 0.5 * median, 2.0 * median

    return lower, upper

def filter_by_length(sequences_dict):
    """Remove obvious length outliers using a robust IQR rule.
    * Records flagged as hypothetical/fragment/partial/putative are discarded first.
    * If <3 sequences remain, no statistical filtering is possible - return them.
    * If ≥5 sequences remain, use Tukey's rule (1.5 x IQR) to define bounds.
    * If 3-4 sequences remain, fall back to a broad window [0.5x, 2x] around the expected length (mode or median).
    """
    # --- step 1: drop low‑confidence annotations -----------------------------
    clean_dict = {}
    for header, seq in sequences_dict.items():
        if not is_suspect_header(header):
            clean_dict[header] = seq
        else:
            print(f"Removed suspect header: {header}")

    if len(clean_dict) < 3:
        print(f"Too few sequences ({len(clean_dict)}) to apply length filter. Returning all non-suspect sequences.")
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
    pass