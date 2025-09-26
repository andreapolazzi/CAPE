from .ncbiprot_api import fetch_protein_ids, fetch_sequences, parse_fasta
import pandas as pd
import regex as re
from typing import List

def is_suspect_header(header: str) -> bool:
    """
    Check if a header contains keywords indicating low confidence annotations.
    
    Args:
        header: FASTA header string
    Returns:
        True if the header is suspect, False otherwise
    """
    lower = header.lower()
    return any(keyword in lower for keyword in ['hypothetical', 'fragment', 'partial', 'putative', 'probable'])   #can be extended

def extract_organism(header: str) -> str:
    """
    Extract the organism name from a FASTA header (mainly from RefSeq entries).
    
    Args:
        header: FASTA header string
    Returns:
        Organism name if found, else None
    """
    match = re.search(r'\[(.*?)\]', header)
    if match:
        return match.group(1)
    return None

def process_ec_number(ec_number: str) -> dict:
        '''
        Process an EC number by fetching sequences and returning a DataFrame and metadata summary.
        
        Args:
            ec_number: Enzyme Commission number as a string
        Returns:
            Dictionary with 'summary' and 'data' (DataFrame)
        '''

        orgs = ['Rhodococcus', 'Nocardiaceae', 'Mycobacteriales', 'Actinomycetes', 'Actinomycetota', 'Bacillati', 'Bacteria']
        for org in orgs:
            protein_ids = fetch_protein_ids(ec_number, org, reviewed=True)
            if protein_ids:
                tax_level = f"{len(protein_ids)} reviewed sequences found at the taxonomy level {org}"
                break
            protein_ids = fetch_protein_ids(ec_number, org, reviewed=False)
            if protein_ids:
                tax_level = f"{len(protein_ids)} non-reviewed sequences found at the taxonomy level {org}"
                break
        else:
            return {
                'summary': None,
                'data': pd.DataFrame(columns=['organism', 'header', 'sequence', 'length'])
            }

        fasta = fetch_sequences(protein_ids)
        parsed_sequences = parse_fasta(fasta)

        hidden_sequences = {}

        # Remove suspect headers
        for header, seq in parsed_sequences.items():
            if not is_suspect_header(header):
                continue
            else:
                hidden_sequences[header] = seq
                del parsed_sequences[header]
                print(f"Hidden suspect entry: {header}")

        # Convert to DataFrame
        records = []
        hidden = []
        for header, sequence in parsed_sequences.items():
            accession_version = re.split(r'\s+', header, maxsplit=1)[0]
            records.append({
                'organism': extract_organism(header),
                'header': header[:header.find('[')] if '[' in header else header,
                'sequence': sequence,
                'length': len(sequence),
                'ncbi_url': f"https://www.ncbi.nlm.nih.gov/protein/{accession_version}"
            })

        for header, sequence in hidden_sequences.items():
            accession_version = re.split(r'\s+', header, maxsplit=1)[0]
            hidden.append({
                'organism': extract_organism(header),
                'header': header[:header.find('[')] if '[' in header else header,
                'sequence': sequence,
                'length': len(sequence),
                'ncbi_url': f"https://www.ncbi.nlm.nih.gov/protein/{accession_version}"
            })

        df = pd.DataFrame(records)
        df_hidden = pd.DataFrame(hidden)
        return {
            'summary': tax_level,
            'data': df,
            'hidden_data': df_hidden if not df_hidden.empty else None
        }

def process_ec_list(ec_list: List[str]) -> List[dict]:
    '''
    Process a list of EC numbers and return results for each.
    
    Args:
        ec_list: List of EC numbers as strings
    Returns:
        List of dictionaries with 'ec_number', 'summary', and 'data' (DataFrame)
    '''
    results = []
    for ec in ec_list:
        result = process_ec_number(ec)
        result['ec_number'] = ec
        results.append(result)
    return results

if __name__ == '__main__':
    test_ecs = ["1.14.12.15", "3.1.1.1"]
    results = process_ec_list(test_ecs)

    print()
    for result in results:
        print(f"EC: {result['ec_number']}")
        print(f"Summary: {result['summary']}")
        print()
        print(result['data'])
        print("-" * 40)
