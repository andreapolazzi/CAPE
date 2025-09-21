#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SelenzymeRF API Client

This script queries the SelenzymeRF API to find enzyme suggestions for chemical reactions
starting from a SMARTS pattern. SelenzymeRF is an updated version of Selenzyme that
uses reaction site mapping to improve enzyme suggestions for biochemical reactions.

Based on:
- SelenzymeRF: updated enzyme suggestion software for unbalanced biochemical reactions
- API endpoint: http://selenzymerf.synbiochem.co.uk
- Original Selenzyme API: http://selenzyme.synbiochem.co.uk/REST

Author: Generated script for SelenzymeRF API usage
"""

import requests
import json
import pandas as pd
import argparse
import sys
from typing import Dict, List, Optional, Any
import time


class SelenzymeRFClient:
    """
    Client for the SelenzymeRF API to query enzyme suggestions from SMARTS patterns.
    """
    
    def __init__(self, base_url: str = "http://selenzymerf.synbiochem.co.uk"):
        """
        Initialize the SelenzymeRF client.
        
        Args:
            base_url: Base URL for the SelenzymeRF API
        """
        self.base_url = base_url.rstrip('/')
        self.rest_endpoint = f"{self.base_url}/REST"
        
    def query_enzymes(self, 
                     smarts: str, 
                     targets: int = 50, 
                     host: Optional[str] = None,
                     noMSA: bool = False,
                     timeout: int = 30) -> Dict[str, Any]:
        """
        Query the SelenzymeRF API for enzyme suggestions based on a SMARTS pattern.
        
        Args:
            smarts: SMARTS pattern representing the chemical reaction
            targets: Maximum number of enzyme results to return (default: 50)
            host: Target host organism for phylogenetic distance calculation (default: '37919' for Rhodococcus opacus)
            noMSA: If True, skip multiple sequence alignment generation
            timeout: Request timeout in seconds
            
        Returns:
            Dictionary containing the API response with enzyme suggestions
            
        Raises:
            requests.RequestException: If the API request fails
            ValueError: If the response format is invalid
        """
        
        # Prepare the request payload
        payload = {
            'smarts': smarts,
            'targets': str(targets),
        }
        
        # Add optional parameters
        if host is not None:
            payload['host'] = host
            
        if noMSA:
            payload['noMSA'] = True
        
        try:
            print(f"Querying SelenzymeRF API with SMARTS: {smarts}")
            
            # Make the API request
            response = requests.post(
                f"{self.rest_endpoint}/Query",
                json=payload,
                timeout=timeout,
                headers={'Content-Type': 'application/json'}
            )
            
            # Check if request was successful
            response.raise_for_status()
            
            # Parse the response
            result = response.json()
            
            if 'data' not in result:
                raise ValueError("Invalid response format: missing 'data' field")
                
            return result
            
        except requests.exceptions.Timeout:
            raise requests.RequestException(f"Request timed out after {timeout} seconds")
        except requests.exceptions.ConnectionError:
            raise requests.RequestException("Failed to connect to SelenzymeRF API")
        except requests.exceptions.HTTPError as e:
            raise requests.RequestException(f"HTTP error {response.status_code}: {e}")
        except json.JSONDecodeError:
            raise ValueError("Invalid JSON response from API")
    
    def parse_enzyme_results(self, api_response: Dict[str, Any]) -> pd.DataFrame:
        """
        Parse the API response into a structured DataFrame.
        
        Args:
            api_response: Raw response from the SelenzymeRF API
            
        Returns:
            pandas.DataFrame with enzyme information
        """
        
        if api_response['data'] is None:
            print("No enzyme data found in response")
            return pd.DataFrame()
        
        try:
            # Parse the data field which contains JSON string
            if isinstance(api_response['data'], str):
                enzyme_data = json.loads(api_response['data'])
            else:
                enzyme_data = api_response['data']
            
            # Check if we have enzyme sequence IDs
            if 'Seq. ID' not in enzyme_data or len(enzyme_data['Seq. ID']) == 0:
                print("No enzyme sequences found")
                return pd.DataFrame()
            
            # Extract enzyme information into rows
            results = []
            seq_ids = enzyme_data['Seq. ID']
            
            # Available columns based on Selenzyme API
            columns = ['Score', 'Seq. ID', 'Description', 'Organism Source', 'Tax. distance', 'EC Number', 'Uniprot protein evidence', 'Consv. Score', 'Combined score']
            
            for seq_id in sorted(seq_ids.keys(), key=lambda x: int(x)):
                row = {'Seq. ID': seq_ids[seq_id]}
                
                # Add other available information
                for col in columns[1:]:  # Skip 'Seq. ID' as it's already added
                    if col in enzyme_data and seq_id in enzyme_data[col]:
                        row[col] = enzyme_data[col][seq_id]
                    else:
                        row[col] = 'N/A'
                
                results.append(row)
            
            return pd.DataFrame(results)
            
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            raise ValueError(f"Error parsing enzyme data: {e}")
    
    def search_enzymes_from_smarts(self, 
                                  smarts: str, 
                                  targets: int = 50,
                                  host: Optional[str] = None,
                                  noMSA: bool = False,
                                  output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Complete workflow: query enzymes and return results as DataFrame.
        
        Args:
            smarts: SMARTS pattern for the reaction
            targets: Maximum number of enzyme results to return
            host: Target host organism for phylogenetic distance calculation
            noMSA: Skip multiple sequence alignment generation
            output_file: Optional file path to save results as CSV
            
        Returns:
            DataFrame with enzyme suggestions
        """
        
        print(f"Starting enzyme search for SMARTS pattern: {smarts}")
        if host:
            print(f"Using host organism: {host}")
        
        # Query the API
        api_response = self.query_enzymes(
            smarts=smarts,
            targets=targets,
            host=host,
            noMSA=noMSA
        )
        
        # Parse results
        results_df = self.parse_enzyme_results(api_response)
        
        if len(results_df) > 0:
            print(f"Found {len(results_df)} enzyme suggestions")
            
            # Save to file if requested
            if output_file:
                results_df.to_csv(output_file, index=False)
                print(f"Results saved to: {output_file}")
        else:
            print("No enzyme suggestions found for this SMARTS pattern")
        
        return results_df


def main():
    """
    Command line interface for the SelenzymeRF API client.
    """
    
    parser = argparse.ArgumentParser(
        description='Query SelenzymeRF API for enzyme suggestions from SMARTS patterns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python selenzymerf_client.py "[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][O:3]"
  python selenzymerf_client.py "CC(=O)O>>CCO" --output results.csv
  python selenzymerf_client.py "c1ccccc1O>>c1ccccc1" --max-results 20
        """
    )
    
    parser.add_argument(
        'smarts',
        help='SMARTS pattern representing the chemical reaction'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Output CSV file to save results'
    )
    
    parser.add_argument(
        '--targets', '-t',
        type=int,
        default=50,
        help='Maximum number of enzyme results to return (default: 50)'
    )
    
    parser.add_argument(
        '--host',
        type=str,
        default='37919',
        help='Target host organism taxonomical ID (default for Rhodococcus opacus: "37919")'
    )
    
    parser.add_argument(
        '--no-msa',
        action='store_true',
        help='Skip multiple sequence alignment generation'
    )
    
    parser.add_argument(
        '--timeout',
        type=int,
        default=30,
        help='Request timeout in seconds (default: 30)'
    )
    
    parser.add_argument(
        '--url',
        default="http://selenzymerf.synbiochem.co.uk",
        help='SelenzymeRF API base URL (default: http://selenzymerf.synbiochem.co.uk)'
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize client
        client = SelenzymeRFClient(base_url=args.url)
        
        # Perform search
        results = client.search_enzymes_from_smarts(
            smarts=args.smarts,
            targets=args.targets,
            host=args.host,
            noMSA=args.no_msa,
            output_file=args.output
        )
        
        # Display results
        if len(results) > 0:
            print("\nTop enzyme suggestions:")
            print("=" * 80)
            print(results.head(10).to_string(index=False))
            
            if len(results) > 10:
                print(f"\n... and {len(results) - 10} more results")
                if args.output:
                    print(f"Full results saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


# Example usage functions
def example_usage():
    """
    Example usage of the SelenzymeRF API client.
    """
    
    # Initialize client
    client = SelenzymeRFClient()
    
    # Example SMARTS patterns for different types of reactions
    examples = [
        # Alcohol oxidation to aldehyde
        "[CH3:1][CH2:2][OH:3]>>[CH3:1][CHO:2]",
        
        # Ester hydrolysis
        "[CH3:1][CO:2][O:3][CH3:4]>>[CH3:1][COOH:2].[OH2:3]",
        
        # Aromatic hydroxylation
        "c1ccccc1>>c1ccc(O)cc1",
        
        # Decarboxylation
        "[CH3:1][COOH:2]>>[CH3:1]"
    ]
    
    for smarts in examples:
        print(f"\n{'='*60}")
        print(f"Searching enzymes for: {smarts}")
        print('='*60)
        
        try:
            results = client.search_enzymes_from_smarts(
                smarts,
                targets=50,
                host="37919"  # Example host organism
            )
            
            if len(results) > 0:
                print(f"Found {len(results)} enzymes:")
                print(results.head().to_string(index=False))
            else:
                print("No enzymes found for this reaction")
                
        except Exception as e:
            print(f"Error processing {smarts}: {e}")
        
        # Small delay between requests to be respectful to the API
        time.sleep(1)


if __name__ == "__main__":
    main()