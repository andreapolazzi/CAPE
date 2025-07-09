import statistics
from ortho_finder import fetch_protein_ids, fetch_sequences, parse_fasta

def infer_expected_length(sequences_dict):
    lengths = [len(seq) for seq in sequences_dict.values()]
    try:
        expected = statistics.mode(lengths)
    except statistics.StatisticsError:
        expected = statistics.median(lengths)
    return expected

def length_score(actual_len, expected_len):
    if actual_len < 0.6 * expected_len:
        return -10  # likely fragment
    elif actual_len > 1.5 * expected_len:
        return -5   # likely fusion
    elif abs(actual_len - expected_len) < 0.15 * expected_len:
        return +5   # ideal
    else:
        return 0
    
def is_suspect_header(header):
    lower = header.lower()
    return any(keyword in lower for keyword in ['hypothetical', 'fragment', 'partial'])   #can be extended

def filter_by_length(sequences_dict):
    expected_len = infer_expected_length(sequences_dict)
    print(f"Inferred expected length: {expected_len}")

    filtered = {}
    for header, sequence in sequences_dict.items():
        if is_suspect_header(header):
            print(f"Removed suspect: {header}")
            continue

        seq_len = len(sequence)
        if 0.6 * expected_len <= seq_len <= 1.5 * expected_len:
            filtered[header] = sequence
        else:
            print(f"Removed outlier: {header[:50]}... | Length: {seq_len}")

    return filtered

if __name__ == '__main__':
    ec = '2.8.3.6'
    org = 'Rhodococcus'
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