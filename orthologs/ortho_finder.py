import requests
from xml.etree import ElementTree

def fetch_protein_ids(ec_number, organism, batch_size=100):
    '''Fetch all protein IDs from NCBI based on EC number and organism'''
    all_ids = []
    start = 0

    # Get total count first
    count_params = {
        'db': 'protein',
        'term': f'{ec_number}[EC/RN Number] AND {organism}[Organism]',
        'retmode': 'xml',
        'retmax': 0
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=count_params)
    tree = ElementTree.fromstring(r.text)
    total_count = int(tree.findtext('.//Count'))
    print(f"Total proteins found for {ec_number} in {organism}: {total_count}")

    # Now fetch in batches
    while start < total_count:
        params = {
            'db': 'protein',
            'term': f'{ec_number}[EC/RN Number] AND {organism}[Organism]',
            'retstart': start,
            'retmax': batch_size
        }
        r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=params)
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
    ec = '1.2.1.10'
    org = 'Rhodococcus'
    protein_ids = fetch_protein_ids(ec, org)
    if not protein_ids:
        print("Nessun ID proteina trovato.")
        exit()
    fasta = fetch_sequences(protein_ids)    #entries as string in FASTA format
    parsed_sequences = parse_fasta(fasta)   #dictionary with headers as keys and sequences as values

    for header, sequence in parsed_sequences.items():
        print(f">{header}\n{sequence}\n")

    print(f"Total sequences fetched: {len(parsed_sequences)}")