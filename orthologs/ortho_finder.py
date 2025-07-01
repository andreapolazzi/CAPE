import requests
from xml.etree import ElementTree

def fetch_protein_ids(ec_number, organism, retmax=100):
    '''Fetch protein IDs from NCBI based on EC number and organism'''
    params = {
        'db': 'protein',
        'term': f'{ec_number}[EC/RN Number] AND {organism}[Organism]',
        'retmax': retmax
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=params)
    tree = ElementTree.fromstring(r.text)
    return [elem.text for elem in tree.findall('.//Id')]

def fetch_sequences(ids):
    '''Fetch protein sequences in FASTA format for a list of IDs'''
    params = {
        'db': 'protein',
        'id': ','.join(ids),
        'rettype': 'fasta',
        'retmode': 'text'
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=params)
    return r.text

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
    ec = '6.2.1.25'
    org = 'Rhodococcus'
    protein_ids = fetch_protein_ids(ec, org)
    if not protein_ids:
        print("Nessun ID proteina trovato.")
        exit()
    fasta = fetch_sequences(protein_ids)    #entries as string in FASTA format
    parsed_sequences = parse_fasta(fasta)   #dictionary with headers as keys and sequences as values

    for header, sequence in parsed_sequences.items():
        print(f">{header}\n{sequence}")