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
from typing import Dict, Optional, Any


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
            noMSA: If True, skip multiple sequence alignment generation (default: False)
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
        else:
            payload['noMSA'] = False
        
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
            api_response: Dictionary containing the API response

        Returns:
            DataFrame with individual enzyme entries as rows and attributes as columns

        Raises:
            ValueError: If the response format is invalid or parsing fails
        
        """
        
        if api_response.get('data') is None:
            print("No enzyme data found in response")
            return pd.DataFrame()
        
        try:
            # Parse the JSON string
            enzyme_data = json.loads(api_response['data'])
            
            if 'Seq. ID' not in enzyme_data:
                print("No sequence IDs found")
                return pd.DataFrame()
            
            # Get all entry indices (they're string numbers)
            entry_indices = set(enzyme_data['Seq. ID'].keys())
            
            # Build rows dynamically
            rows = []
            for idx in sorted(entry_indices, key=int):  # Sort numerically
                row = {}
                
                # Extract data for each column
                for col_name, col_data in enzyme_data.items():
                    if idx in col_data and col_data[idx] is not None:
                        row[col_name] = col_data[idx]
                    else:
                        row[col_name] = None  # Keep nulls as None instead of 'N/A'
                
                rows.append(row)
            
            return pd.DataFrame(rows)
            
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            raise ValueError(f"Error parsing enzyme data: {e}")
    
    
    def search_enzymes(self, smarts: str, **kwargs) -> pd.DataFrame:
        """
        Main method to search for enzymes
        
        Args:
            smarts: SMARTS pattern representing the chemical reaction
            targets: Maximum number of enzyme results to return (default: 50)
            host: Target host organism for phylogenetic distance calculation (default: '37919' for Rhodococcus opacus)
            noMSA: If True, skip multiple sequence alignment generation (default: False)
            timeout: Request timeout in seconds

        Returns:
            DataFrame with individual enzyme entries as rows and attributes as columns

        Raises:
            requests.RequestException: If the API request fails
            ValueError: If the response format is invalid or parsing fails
        """
        
        api_response = self.query_enzymes(smarts, **kwargs)
        df = self.parse_enzyme_results(api_response)
        
        return df


if __name__ == "__main__":
    pass