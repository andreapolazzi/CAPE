"""
Batch SelenzymeRF API Client

Simple batch processor for SelenzymeRF API that takes a list of SMARTS strings
and returns a dictionary mapping each SMARTS to its results DataFrame.

"""

import time
from typing import List, Dict, Optional
import pandas as pd
from .api_selenzyme import SelenzymeRFClient

def process_smarts(
    smarts: str,
    targets: int = 50,
    host: Optional[str] = "37919",
    api_url: str = "http://selenzymerf.synbiochem.co.uk",
    timeout: int = 30
) -> pd.DataFrame:
    """
    Process a single SMARTS string using the SelenzymeRF API.
    
    Args:
        smarts: SMARTS pattern to process
        targets: Maximum number of enzyme results (default: 50)
        host: Host organism taxonomic ID (default: "37919" for Rhodococcus opacus)
        api_url: SelenzymeRF API base URL
        timeout: Request timeout in seconds
        
    Returns:
        DataFrame containing enzyme suggestions
    """
    
    # Initialize the API client
    client = SelenzymeRFClient(base_url=api_url)
    
    # Query the API
    df = client.search_enzymes(
        smarts=smarts,
        targets=targets,
        host=host,
        noMSA=True,
        timeout=timeout
    )
    
    return df


def process_smarts_batch(
    smarts_list: List[str],
    targets: int = 50,
    host: Optional[str] = "37919",
    delay_between_requests: float = 0.5,
    api_url: str = "http://selenzymerf.synbiochem.co.uk",
    timeout: int = 60
) -> Dict[str, pd.DataFrame]:
    """
    Process a batch of SMARTS strings using the SelenzymeRF API.
    
    Args:
        smarts_list: List of SMARTS patterns to process
        targets: Maximum number of enzyme results per query
        host: Host organism taxonomic ID (default: "37919" for Rhodococcus opacus)
        delay_between_requests: Seconds to wait between API calls
        api_url: SelenzymeRF API base URL
        timeout: Request timeout in seconds
        
    Returns:
        Dictionary mapping each SMARTS string to its results DataFrame
    """
    
    # Initialize the API client
    client = SelenzymeRFClient(base_url=api_url)
    
    results = {}
    
    for i, smarts in enumerate(smarts_list):
        print(f"Processing SMARTS {i+1}/{len(smarts_list)}: {smarts[:50]}...")
        
        try:
            # Query the API
            df = client.search_enzymes(
                smarts=smarts,
                targets=targets,
                host=host,
                noMSA=True,
                timeout=timeout
            )
            
            results[smarts] = df
            print(f"  Found {len(df)} enzyme suggestions")
            
        except Exception as e:
            print(f"  Error: {e}")
            # Store empty DataFrame on error
            results[smarts] = pd.DataFrame()
        
        # Rate limiting - wait between requests except for the last one
        if i < len(smarts_list) - 1:
            time.sleep(delay_between_requests)
    
    return results


def process_smarts_filtered(
    smarts: str,
    include_columns: Optional[List[str]] = None,
    **kwargs
) -> pd.DataFrame:
    """
    Process a single SMARTS string and return only specified columns.
    Args:
        smarts: SMARTS pattern to process
        include_columns: List of column names to include in results
        **kwargs: Additional arguments passed to process_smarts
    Returns:
        DataFrame containing filtered enzyme suggestions
    """
    # Get full results
    df = process_smarts(smarts, **kwargs)
    
    # Filter columns if specified
    if include_columns is not None and len(df) > 0:
        # Only include columns that exist in the DataFrame
        available_columns = [col for col in include_columns if col in df.columns]
        if available_columns:
            return df[available_columns]
        else:
            return pd.DataFrame()
    
    return df

def process_smarts_batch_filtered(
    smarts_list: List[str],
    include_columns: Optional[List[str]] = None,
    **kwargs
) -> Dict[str, pd.DataFrame]:
    """
    Process a batch of SMARTS strings and return only specified columns.
    
    Args:
        smarts_list: List of SMARTS patterns to process
        include_columns: List of column names to include in results
        **kwargs: Additional arguments passed to process_smarts_batch
        
    Returns:
        Dictionary mapping each SMARTS string to its filtered results DataFrame
    """
    
    # Get full results
    results = process_smarts_batch(smarts_list, **kwargs)
    
    # Filter columns if specified
    if include_columns is not None:
        filtered_results = {}
        for smarts, df in results.items():
            if len(df) > 0:
                # Only include columns that exist in the DataFrame
                available_columns = [col for col in include_columns if col in df.columns]
                if available_columns:
                    filtered_results[smarts] = df[available_columns]
                else:
                    filtered_results[smarts] = pd.DataFrame()
            else:
                filtered_results[smarts] = df
        return filtered_results
    
    return results


if __name__ == "__main__":
    # Example usage for batch processing
    # test_smarts = [
    #     "CC(=O)O>>CCO",  # Simple ester reduction
    #     "c1ccccc1O>>c1ccc(O)cc1",  # Phenol hydroxylation
    #     "[H]C1=CC(=CC=C1C=O)O[H]>>[H]C1=CC(=CC=C1C(=O)O[H])O[H]"  # Aldehyde oxidation
    # ]
    
    # print("Testing batch processing...")
    # results = process_smarts_batch(test_smarts, targets=10)
    
    # print("\nResults summary:")
    # for smarts, df in results.items():
    #     print(f"{smarts}: {len(df)} enzymes found")
    #     if len(df) > 0:
    #         print(f"  Top enzyme: {df.iloc[0]['Description']}")
    
    # print("\nTesting with column filtering...")
    # filtered_results = process_smarts_batch_filtered(
    #     test_smarts[:2],  # Just first two for demo
    #     include_columns=["Seq. ID", "Description", "Organism Source", "EC Number", "Score"],
    #     targets=5
    # )
    
    # for smarts, df in filtered_results.items():
    #     print(f"\nFiltered results for: {smarts}")
    #     if len(df) > 0:
    #         print(df.head())
    #     else:
    #         print("No results")

    # Example usage for single SMARTS processing
    test_smarts_single = "CC(=O)O>>CCO"  # Simple ester reduction
    print("Testing single SMARTS processing...")
    df_single = process_smarts(test_smarts_single, targets=10)
    print(df_single.head())