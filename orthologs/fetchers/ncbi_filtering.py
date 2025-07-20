import requests
from xml.etree import ElementTree
import sys

def fetch_protein_ids(ec_number, organism, batch_size=100, reviewed=False):
    '''Fetch all protein IDs from NCBI based on EC number and organism, optionally filtering for Swiss-Prot'''
    all_ids = []
    start = 0

    base_query = f'{ec_number}[ECNO] AND {organism}[ORGN]'
    if reviewed:
        base_query += ' AND srcdb_swiss-prot[PROP]'

    # Get total count first
    count_params = {
        'db': 'protein',
        'term': base_query,
        'retmode': 'xml',
        'retmax': 0
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=count_params)
    r.raise_for_status()
    tree = ElementTree.fromstring(r.text)
    total_count = int(tree.findtext('.//Count'))
    if reviewed:
        print(f"Total reviwed proteins found for {ec_number} in {organism}: {total_count}")
    else:
        print(f"Total non-reviwed proteins found for {ec_number} in {organism}: {total_count}")

    # Now fetch in batches
    while start < total_count:
        params = {
            'db': 'protein',
            'term': base_query,
            'retstart': start,
            'retmax': batch_size
        }
        r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=params)
        r.raise_for_status()
        tree = ElementTree.fromstring(r.text)
        ids = [elem.text for elem in tree.findall('.//Id')]
        all_ids.extend(ids)
        start += batch_size

    return all_ids

def fetch_sequences(ids, chunk_size=200):
    '''Fetch protein sequences in FASTA format for a list of IDs, in chunks'''
    fasta_all = ''
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        params = {
            'db': 'protein',
            'id': ','.join(chunk),
            'rettype': 'fasta',
            'retmode': 'text'
        }
        r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=params)
        r.raise_for_status()
        fasta_all += r.text
    return fasta_all

def parse_fasta(fasta_text):
    '''Parse FASTA formatted text into a dictionary'''
    sequences = {}
    header = None
    seq_lines = []
    for line in fasta_text.splitlines():
        if line.startswith('>'):
            if header is not None:
                sequences[header] = ''.join(seq_lines)
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        sequences[header] = ''.join(seq_lines)
    return sequences

if __name__ == '__main__':
    from remove_outliers import filter_by_length
    ec = '3.7.1.21'
    orgs = ['Rhodococcus', 'Nocardiaceae', 'Mycobacteriales', 'Actinomycetes', 'Actinomycetota', 'Bacillati']

    for org in orgs:
        print(f"\nFetching protein IDs for EC {ec} in {org}...")
        protein_ids = fetch_protein_ids(ec, org, reviewed=True)
        if protein_ids:
            print(f"Found {len(protein_ids)} reviewed protein IDs.")
            break

        print("No reviewed protein IDs found, trying without reviewed filter...")
        protein_ids = fetch_protein_ids(ec, org, reviewed=False)
        if protein_ids:
            print(f"Found {len(protein_ids)} non‑reviewed protein IDs.")
            break

        print("No protein IDs found for this taxon.")
    else:
        print("No protein IDs found in any taxonomic level.")
        sys.exit(1)

    fasta = fetch_sequences(protein_ids)    #entries as string in FASTA format
    parsed_sequences = parse_fasta(fasta)   #dictionary with headers as keys and sequences as values

    filtered_sequences = filter_by_length(parsed_sequences)
    print(f"\nRemaining sequences after outlier removal: {len(filtered_sequences)}")
    for header, sequence in filtered_sequences.items():
        print(f">{header}\n{sequence}\n")

    print(f"Total sequences fetched: {len(parsed_sequences)}; kept after filtering: {len(filtered_sequences)}")