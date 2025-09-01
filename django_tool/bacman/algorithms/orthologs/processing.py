from filters.ncbi_filtering import fetch_protein_ids, fetch_sequences, parse_fasta
import pandas as pd
import regex as re

def is_suspect_header(header):
    """Check if a header contains keywords indicating low confidence annotations."""
    lower = header.lower()
    return any(keyword in lower for keyword in ['hypothetical', 'fragment', 'partial', 'putative', 'probable'])   #can be extended

def extract_organism(header):
    match = re.search(r'\[(.*?)\]', header)
    if match:
        return match.group(1)
    return None

def process_ec_number(ec_number):
        '''Process an EC number by fetching sequences and returning a DataFrame and metadata summary'''

        orgs = ['Rhodococcus', 'Nocardiaceae', 'Mycobacteriales', 'Actinomycetes', 'Actinomycetota', 'Bacillati']
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
        for header, seq in parsed_sequences.items():
            if not is_suspect_header(header):
                parsed_sequences[header] = seq
            else:
                print(f"Removed suspect header: {header}")

        records = []
        for header, sequence in parsed_sequences.items():
            records.append({
                'organism': extract_organism(header),
                'header': header[:header.find('[')] if '[' in header else header,
                'sequence': sequence,
                'length': len(sequence)
            })

        df = pd.DataFrame(records)
        return {
            'summary': tax_level,
            'data': df
        }

def process_ec_list(ec_list):
    '''Process a list of EC numbers and return results for each'''
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
