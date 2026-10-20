import requests
from xml.etree import ElementTree
from typing import List
import time

def fetch_protein_ids(ec_number: str, organism: str, batch_size=100, reviewed=False) -> list:
    '''
    Fetch all protein IDs from NCBI based on EC number and organism, optionally filtering for Swiss-Prot

    Args:
        ec_number: Enzyme Commission number associated with the protein
        organism: Organism taxonomic ID (as found in https://www.ncbi.nlm.nih.gov/taxonomy)
        batch_size: Number of IDs to fetch per request (default: 100)
        reviewed: If True, only fetch reviewed (Swiss-Prot) entries (default: False)
    Returns:
        List of protein IDs
    '''
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
    time.sleep(0.34)  # To respect NCBI's rate limit of 3 requests per second
    r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=count_params)
    r.raise_for_status()
    tree = ElementTree.fromstring(r.text)
    total_count = int(tree.findtext('.//Count'))
    if reviewed:
        print(f"Total reviewed proteins found for {ec_number} in {organism}: {total_count}")
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
        time.sleep(0.34)  # To respect NCBI's rate limit of 3 requests per second
        r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=params)
        r.raise_for_status()
        tree = ElementTree.fromstring(r.text)
        ids = [elem.text for elem in tree.findall('.//Id')]
        all_ids.extend(ids)
        start += batch_size

    return all_ids

def fetch_sequences(ids: List[str], chunk_size: int = 150) -> str:
    '''
    Fetch protein sequences in FASTA format for a list of IDs, in chunks
    
    Args:
        ids: List of protein IDs
        chunk_size: Number of IDs to fetch per request (default: 150)
    Returns:
        Concatenated FASTA formatted string of all sequences
    '''
    fasta_all = ''
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i:i+chunk_size]
        params = {
            'db': 'protein',
            'id': ','.join(chunk),
            'rettype': 'fasta',
            'retmode': 'text'
        }
        time.sleep(0.34)  # To respect NCBI's rate limit of 3 requests per second
        r = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi', params=params)
        r.raise_for_status()
        fasta_all += r.text
    return fasta_all

def parse_fasta(fasta_text: str) -> dict:
    '''
    Parse FASTA formatted text into a dictionary
    
    Args:
        fasta_text: FASTA formatted string
    Returns:
        Dictionary mapping headers to sequences
    '''
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
    pass