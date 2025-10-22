try:
    from .ncbiprot_api import fetch_protein_ids, fetch_sequences, parse_fasta
except ImportError:
    from ncbiprot_api import fetch_protein_ids, fetch_sequences, parse_fasta
import pandas as pd
import regex as re
from typing import List, Union
from pathlib import Path
import sys
import argparse

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
            Dictionary with 'summary', 'data' (DataFrame) and 'hidden_data' (DataFrame of suspect entries)
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
        suspect_headers = []
        # Remove suspect headers
        for header, seq in parsed_sequences.items():
            if is_suspect_header(header):
                hidden_sequences[header] = seq
                suspect_headers.append(header)
                print(f"Hidden suspect entry: {header}")

        for header in suspect_headers:
            del parsed_sequences[header]

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
        List of dictionaries with 'ec_number', 'summary', 'data' (DataFrame) and 'hidden_data' (DataFrame of suspect entries)
    '''
    results = []
    for ec in ec_list:
        result = process_ec_number(ec)
        result['ec_number'] = ec
        results.append(result)
    return results

def save_results_to_csv(results: List[dict], output_dir: Union[Path, str], include_hidden: bool = False):
    """Save results to CSV files."""
    output_dir.mkdir(exist_ok=True)
    
    for result in results:
        ec_number = result['ec_number']
        safe_ec = ec_number.replace('.', '_')
        
        # Save main data
        if not result['data'].empty:
            main_file = output_dir / f"EC_{safe_ec}_sequences.csv"
            result['data'].to_csv(main_file, index=False)
            print(f"Saved main sequences to: {main_file}")
        
        # Save hidden data if requested and exists
        if include_hidden and result.get('hidden_data') is not None:
            hidden_file = output_dir / f"EC_{safe_ec}_hidden.csv"
            result['hidden_data'].to_csv(hidden_file, index=False)
            print(f"Saved hidden sequences to: {hidden_file}")

def read_ec_list_from_file(file_path: Path) -> List[str]:
    """Read EC numbers from a file (one per line)."""
    try:
        with open(file_path, 'r') as f:
            ec_numbers = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        return ec_numbers
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Process EC numbers and fetch protein sequences from NCBI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 1.14.12.15                           # Process single EC number
  %(prog)s 1.14.12.15 3.1.1.1                  # Process multiple EC numbers
  %(prog)s --file ec_numbers.txt                # Process EC numbers from file
  %(prog)s 1.14.12.15 --outdir results/        # Save to CSV files
  %(prog)s 1.14.12.15 --include-hidden         # Include suspect sequences in output
        """
    )
    
    # Input options
    parser.add_argument(
        'ec_numbers',
        nargs='*',
        help='EC numbers to process (e.g., 1.14.12.15 3.1.1.1)'
    )
    parser.add_argument(
        '--file', '-f',
        type=Path,
        help='File containing EC numbers (one per line, # for comments)'
    )
    
    # Output options
    parser.add_argument(
        '--outdir', '-o',
        type=Path,
        help='Output directory for CSV files (default: print to console)'
    )

    parser.add_argument(
        '--include-hidden',
        action='store_true',
        help='Include suspect sequences (hypothetical, fragment, etc.) in output files'
    )
    
    args = parser.parse_args()
    
    # Get EC numbers
    if args.file:
        ec_numbers = read_ec_list_from_file(args.file)
    else:
        ec_numbers = args.ec_numbers
    
    if not ec_numbers:
        print("No EC numbers provided.")
        sys.exit(1)
    
    try:
        results = process_ec_list(ec_numbers)
    except Exception as e:
        print(f"Error processing EC numbers: {e}")
        sys.exit(1)
    
    # Save results
    if args.outdir:
        save_results_to_csv(results, args.outdir, args.include_hidden)

    # Print results to console if no output directory specified
    else:
        for result in results:
            ec_number = result['ec_number']
            print(f"\nEC Number: {ec_number}")
            if result['summary']:
                print(f"Summary: {result['summary']}")
            if not result['data'].empty:
                print("Main Sequences:")
                print(result['data'].head(min(5, len(result['data']))))
            else:
                print("No sequences found.")
            if result.get('hidden_data') is not None:
                print("\nHidden Sequences:")
                print(result['hidden_data'].head(min(5, len(result['hidden_data']))))


if __name__ == '__main__':
    main()