import statistics
from fetchers.ncbi_filtering import fetch_protein_ids, fetch_sequences, parse_fasta

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
    if len(sequences_dict) < 3:
        print("Too few sequences to apply length filter. Returning all non-suspect sequences.")
        return {
            header: sequence
            for header, sequence in sequences_dict.items()
            if not is_suspect_header(header)
        }

    expected_len = infer_expected_length(sequences_dict)
    print(f"Inferred expected length: {expected_len}")

    filtered = {}
    relaxed = len(sequences_dict) < 10
    if relaxed:
        print("Applying relaxed filtering due to small number of sequences.")

    for header, sequence in sequences_dict.items():
        if is_suspect_header(header):
            print(f"Removed suspect: {header}")
            continue

        seq_len = len(sequence)
        if relaxed:
            if 0.5 * expected_len <= seq_len <= 2.0 * expected_len:
                filtered[header] = sequence
            else:
                print(f"Removed extreme outlier: {header[:50]}... | Length: {seq_len}")
        else:
            if 0.6 * expected_len <= seq_len <= 1.5 * expected_len:
                filtered[header] = sequence
            else:
                print(f"Removed outlier: {header[:50]}... | Length: {seq_len}")

    return filtered

if __name__ == '__main__':
    ec = '1.13.11.8'
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