import requests
from xml.etree import ElementTree

def fetch_protein_ids(ec_number, organism, retmax=100):
    params = {
        'db': 'protein',
        'term': f'{ec_number}[EC] AND {organism}[Organism]',
        'retmax': retmax
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=params)
    tree = ElementTree.fromstring(r.text)
    return [elem.text for elem in tree.findall('.//Id')]

def fetch_sequences(ids):
    params = {
        'db': 'protein',
        'id': ','.join(ids),
        'rettype': 'fasta',
        'retmode': 'text'
    }
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=params)
    return r.text

if __name__ == '__main__':
    ec = '2.8.3.6'
    org = 'Rhodococcus'
    protein_ids = fetch_protein_ids(ec, org)
    if not protein_ids:
        print("Nessun ID proteina trovato.")
        exit()
    fasta = fetch_sequences(protein_ids)
    print(fasta)