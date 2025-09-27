#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NOEC SelenzymeRF API Interface

Main interface for SelenzymeRF API operations. Can be used both as a standalone
CLI tool and imported as a module for backend integration.

This module provides:
- Single SMARTS processing
- Batch SMARTS processing with rate limiting
- Command-line interface
- Backend-friendly functions for integration

"""

import argparse
import sys
import time
from typing import List, Dict, Optional, Union
import pandas as pd
from pathlib import Path

from .api_selenzyme import SelenzymeRFClient

def query_smarts(
    smarts: Union[str, List[str]],
    targets: int = 50,
    host: Optional[str] = "37919",
    delay_between_requests: float = 1,
    timeout: int = 60,
    api_url: str = "http://selenzymerf.synbiochem.co.uk"
) -> Dict[str, pd.DataFrame]:
    """
    Query SelenzymeRF API for one or more SMARTS patterns.
    
    Unified function that handles both single and batch processing automatically.
    Always continues on errors and returns what it can.
    
    Args:
        smarts: Single SMARTS string or list of SMARTS patterns
        targets: Maximum number of enzyme results per query (default: 50)
        host: Host organism taxonomic ID (default: "37919" for Rhodococcus opacus)
        delay_between_requests: Seconds to wait between batch requests (default: 1)
        timeout: Request timeout in seconds (default: 60)
        api_url: SelenzymeRF API base URL
        
    Returns:
        Dictionary mapping each SMARTS string to its results DataFrame
        For single SMARTS input, returns dict with one key-value pair
        Failed queries return empty DataFrames
    """
    # Normalize input to list
    if isinstance(smarts, str):
        smarts_list = [smarts]
        is_single = True
    else:
        smarts_list = smarts
        is_single = False
    
    client = SelenzymeRFClient(base_url=api_url)
    results = {}
    
    for i, smarts_pattern in enumerate(smarts_list):
        try:
            df = client.search_enzymes(
                smarts=smarts_pattern,
                targets=targets,
                host=host,
                noMSA=True,
                timeout=timeout
            )
            results[smarts_pattern] = df
            
        except Exception as e:
            # On error, store empty DataFrame
            results[smarts_pattern] = pd.DataFrame()
            print(f"Warning: Failed to process SMARTS '{smarts_pattern}': {e}")
        
        # Rate limiting - wait between requests except for the last one or single queries
        if not is_single and i < len(smarts_list) - 1:
            time.sleep(delay_between_requests)
    
    return results


def get_single_result(results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Helper function to extract single DataFrame from query_smarts results.
    
    Args:
        results: Results dictionary from query_smarts
        
    Returns:
        First DataFrame from results
        
    Raises:
        ValueError: If results is empty
    """
    if not results:
        raise ValueError("No results to extract")
    
    return next(iter(results.values()))


def save_results_to_files(
    results: Dict[str, pd.DataFrame],
    outdir: Union[str, Path],
    prefix: str = "selenzymeRF_"
) -> List[str]:
    """
    Save batch results to separate CSV files.
    
    Args:
        results: Dictionary from query_batch_smarts
        outdir: Directory to save files
        prefix: Filename prefix
        
    Returns:
        List of created file paths
    """
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)
    
    saved_files = []
    
    for i, (smarts, df) in enumerate(results.items()):
        if len(df) > 0:
            filename = f"{prefix}{i+1:03d}.csv"
            filepath = outdir / filename
            df.to_csv(filepath, index=False)
            saved_files.append(str(filepath))
    
    return saved_files


def load_smarts_from_file(filepath: str) -> List[str]:
    """Load SMARTS patterns from a text file (one per line)."""
    with open(filepath, 'r') as f:
        smarts_list = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    return smarts_list


def main():
    """Command line interface for SelenzymeRF operations."""
    
    parser = argparse.ArgumentParser(
        description='NOEC SelenzymeRF API Interface - Query enzymes from SMARTS patterns',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single SMARTS query
  python noec_api.py "CC(=O)O>>CCO" --output results.csv
  
  # Multiple SMARTS as arguments
  python noec_api.py "CC(=O)O>>CCO" "c1ccccc1O>>c1ccccc1" --output-dir ./results/
  
  # From file
  python noec_api.py --file smarts_list.txt --outdir ./results/
        """
    )
    
    # Input options
    parser.add_argument(
        'smarts',
        nargs='*',
        help='SMARTS pattern(s) or path to file containing SMARTS (one per line)'
    )
    parser.add_argument(
        '--file', '-f',
        help='Path to file containing SMARTS patterns (one per line)'
    )
    
    # Output options
    parser.add_argument(
        '--output', '-o',
        help='Output CSV file for single query'
    )
    parser.add_argument(
        '--outdir',
        help='Output directory for multiple results (creates separate files)'
    )
    
    # API parameters
    parser.add_argument(
        '--targets', '-t',
        type=int,
        default=50,
        help='Maximum number of enzyme results per query (default: 50)'
    )
    parser.add_argument(
        '--host',
        default="37919",
        help='Host organism taxonomic ID (default: "37919" for Rhodococcus opacus)'
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=30,
        help='Request timeout in seconds (default: 30)'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=1,
        help='Delay between requests when processing multiple SMARTS (default: 1)'
    )
    parser.add_argument(
        '--api-url',
        default="http://selenzymerf.synbiochem.co.uk",
        help='SelenzymeRF API base URL'
    )
    
    args = parser.parse_args()
    
    try:
        # Determine input source and get SMARTS list
        smarts_list = []
        
        # From file
        if args.file:
            smarts_list = load_smarts_from_file(args.file)
            print(f"Loaded {len(smarts_list)} SMARTS patterns from {args.file}")
        
        # From command line arguments
        elif args.smarts:
            # Check if it's a single file path
            if len(args.smarts) == 1 and Path(args.smarts[0]).exists():
                smarts_list = load_smarts_from_file(args.smarts[0])
                print(f"Loaded {len(smarts_list)} SMARTS patterns from {args.smarts[0]}")
            else:
                # Direct SMARTS patterns
                smarts_list = args.smarts
        
        # No input provided
        else:
            parser.error("No SMARTS pattern(s) provided. Use positional arguments or --file.")
        
        print(f"Processing {len(smarts_list)} SMARTS pattern(s)...")
        
        results = query_smarts(
            smarts=smarts_list,
            targets=args.targets,
            host=args.host,
            delay_between_requests=args.delay,
            timeout=args.timeout,
            api_url=args.api_url
        )
        
        # Handle output
        if args.output or args.outdir:
            if args.output and len(smarts_list) > 1:
                print("Warning: --output specified but multiple SMARTS provided. Using --outdir instead.")
                args.output = None
            
            if args.output:
                # Single output file
                df = get_single_result(results)
                if len(df) > 0:
                    df.to_csv(args.output, index=False)
                    print(f"Results saved to {args.output}")
                else:
                    print(f"No results to save for the provided SMARTS.")
            
            elif args.outdir:
                # Multiple output files
                saved_files = save_results_to_files(results, args.outdir)
                if saved_files:
                    print(f"Results saved to {len(saved_files)} files in {args.outdir}")
                else:
                    print("No results to save for the provided SMARTS.")
        else:
            # Print summary to console
            for smarts, df in results.items():
                print(f"\nSMARTS: {smarts}")
                if len(df) > 0:
                    print(df.head())
                    print(f"... ({len(df)} total results)")
                else:
                    print("No results found.")
            
    
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()