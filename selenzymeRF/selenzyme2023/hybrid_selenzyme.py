#!/usr/bin/env python3
"""
Hybrid SelenzymeRF Client - API with Local Fallback

This script provides a robust SelenzymeRF interface that:
1. First tries the SelenzymeRF API (fast, no local setup needed)
2. Falls back to local execution if API fails (reliable, full control)

Requirements:
- pandas
- requests (for API)
- Local SelenzymeRF installation (for fallback)

Example:
    client = HybridSelenzymeRF(
        api_url="http://selenzymerf.synbiochem.co.uk",
        repo_root="/work/igem/selenzyme2023"
    )
    
    df = client.search_enzymes(
        smarts="CC(=O)O>>CCO",
        targets=50,
        host_taxid=37919
    )
"""

from __future__ import annotations
import os, sys, subprocess
from pathlib import Path
import pandas as pd
import requests
import json
import time
from typing import Optional, Dict, Union
import warnings


class HybridSelenzymeRF:
    """
    Hybrid SelenzymeRF client that tries API first, then falls back to local execution.
    """
    
    def __init__(self, 
                 api_url: str = "http://selenzymerf.synbiochem.co.uk",
                 repo_root: Optional[str] = None,
                 datadir: Optional[str] = None,
                 outdir: Optional[str] = None,
                 python_exe: Optional[str] = None,
                 api_timeout: int = 30,
                 local_timeout: Optional[int] = None):
        """
        Initialize the hybrid client.
        
        Args:
            api_url: Base URL for SelenzymeRF API
            repo_root: Path to local SelenzymeRF repository
            datadir: Data directory for local execution
            outdir: Output directory for local execution  
            python_exe: Python executable for local execution
            api_timeout: Timeout for API requests (seconds)
            local_timeout: Timeout for local execution (seconds)
        """
        self.api_url = api_url.rstrip('/') if api_url else None
        self.repo_root = Path(repo_root).resolve() if repo_root else None
        self.datadir = datadir
        self.outdir = outdir
        self.python_exe = python_exe
        self.api_timeout = api_timeout
        self.local_timeout = local_timeout
        
        # Validate local setup if repo_root provided
        if self.repo_root:
            self._validate_local_setup()
    
    def _validate_local_setup(self):
        """Validate that local SelenzymeRF is properly set up."""
        selenzy_path = self.repo_root / "gitcode2023" / "selenzyPro" / "selenzy.py"
        if not selenzy_path.is_file():
            warnings.warn(f"Local selenzy.py not found at {selenzy_path}. Local fallback disabled.")
            self.repo_root = None
    
    def search_enzymes(self,
                      smarts: str,
                      targets: int = 50,
                      host_taxid: Optional[int] = None,
                      nomsa: bool = True,
                      outfile: str = "selenzyme_results.csv") -> pd.DataFrame:
        """
        Search for enzymes using hybrid approach (API first, then local fallback).
        
        Args:
            smarts: SMARTS pattern for the reaction
            targets: Maximum number of enzyme results
            host_taxid: Host organism taxonomic ID (e.g., 37919 for Rhodococcus opacus)
            nomsa: Skip multiple sequence alignment
            outfile: Output filename for local execution
            
        Returns:
            DataFrame with enzyme suggestions
        """
        
        # Method 1: Try API first
        # if self.api_url:
        #     try:
        #         print("🌐 Trying SelenzymeRF API...")
        #         result = self._search_via_api(smarts, targets, host_taxid, nomsa)
        #         print("✅ API query successful!")
        #         return result
        #     except Exception as e:
        #         print(f"❌ API failed: {e}")
        #         if self.repo_root is None:
        #             raise RuntimeError(f"API failed and no local fallback available: {e}")
        #         print("🔄 Falling back to local execution...")
        
        # Method 2: Local execution fallback
        if self.repo_root:
            try:
                print("💻 Running local SelenzymeRF...")
                result = self._search_locally(smarts, targets, host_taxid, outfile)
                print("✅ Local execution successful!")
                return result
            except Exception as e:
                print(f"❌ Local execution failed: {e}")
                raise RuntimeError(f"Both API and local execution failed. Local error: {e}")
        
        # If we get here, no methods are available
        raise RuntimeError("No execution methods available (both API URL and repo_root missing)")
    
    def _search_via_api(self, smarts: str, targets: int, host_taxid: Optional[int], nomsa: bool) -> pd.DataFrame:
        """Execute query via SelenzymeRF API."""
        
        # Prepare payload
        payload = {
            'smarts': smarts,
            'targets': str(targets),
        }
        
        if host_taxid is not None:
            payload['host'] = str(host_taxid)
        
        if nomsa:
            payload['noMSA'] = True
        
        # Make API request
        response = requests.post(
            f"{self.api_url}/REST/Query",
            json=payload,
            timeout=self.api_timeout,
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        
        # Parse response
        result = response.json()
        if 'data' not in result:
            raise ValueError("Invalid API response: missing 'data' field")
        
        if result['data'] is None:
            print("No enzyme data found in API response")
            return pd.DataFrame()
        
        # Parse enzyme data
        if isinstance(result['data'], str):
            enzyme_data = json.loads(result['data'])
        else:
            enzyme_data = result['data']
        
        if 'Seq. ID' not in enzyme_data or len(enzyme_data['Seq. ID']) == 0:
            print("No enzyme sequences found")
            return pd.DataFrame()
        
        # Convert to DataFrame
        results = []
        seq_ids = enzyme_data['Seq. ID']
        
        # Standard columns from API
        possible_columns = [
            'Score', 'Description', 'Organism Source', 'Tax. distance', 
            'EC Number', 'Uniprot protein evidence', 'Consv. Score', 'Combined score'
        ]
        
        for seq_id in sorted(seq_ids.keys(), key=lambda x: int(x)):
            row = {'Seq. ID': seq_ids[seq_id]}
            
            for col in possible_columns:
                if col in enzyme_data and seq_id in enzyme_data[col]:
                    row[col] = enzyme_data[col][seq_id]
                else:
                    row[col] = 'N/A'
            
            results.append(row)
        
        df = pd.DataFrame(results)
        df['source'] = 'API'
        return df
    
    def _search_locally(self, smarts: str, targets: int, host_taxid: Optional[int], 
                       outfile: str) -> pd.DataFrame:
        """Execute query via local SelenzymeRF installation."""
        
        # Locate selenzy.py
        selenzy_path = self.repo_root / "gitcode2023" / "selenzyPro" / "selenzy.py"
        if not selenzy_path.is_file():
            raise FileNotFoundError(f"selenzy.py not found at: {selenzy_path}")
        
        # Set up directories
        if self.datadir is None:
            datadir = self.repo_root / "gitcode2023" / "selenzyPro" / "data" / "data_2023"
        else:
            datadir = Path(self.datadir)
        
        if self.outdir is None:
            outdir = self.repo_root / "runout"
        else:
            outdir = Path(self.outdir)
        
        # Create output directory
        outdir.mkdir(parents=True, exist_ok=True)
        outdir_with_sep = str(outdir) + os.sep
        
        # Build command as required by selenzy.py
        py = self.python_exe or sys.executable
        cmd = [py, str(selenzy_path),
               "-smarts", smarts,
               str(datadir),
               outdir_with_sep,
               "-tar", str(targets)]
        
        # Avoid MSA to save time when running locally
        cmd += ["-NoMSA"]
        
        if host_taxid is not None:
            cmd += ["-host", str(host_taxid)]
        cmd += ["-outfile", outfile]
        
        # Execute
        proc = subprocess.run(
            cmd, cwd=selenzy_path.parent, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=self.local_timeout
        )
        
        if proc.returncode != 0:
            raise RuntimeError(f"Local SelenzymeRF failed (code {proc.returncode}):\n"
                              f"CMD: {' '.join(cmd)}\n"
                              f"STDERR: {proc.stderr}")
        
        # Find and load CSV
        csv_path = outdir / outfile

        if not csv_path.exists():
            raise FileNotFoundError(f"Expected CSV not found: {outdir/outfile}")
        
        df = pd.read_csv(csv_path)
        df['source'] = 'Local'

        return df



def main():
    """Command line interface for hybrid SelenzymeRF."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Hybrid SelenzymeRF client (API + Local fallback)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Try API first, fallback to local
  python hybrid_selenzyme.py "CC(=O)O>>CCO" --repo-root /path/to/selenzyme
  
  # API only (no local repo specified)
  python hybrid_selenzyme.py "CC(=O)O>>CCO"
  
  # With specific host organism and targets
  python hybrid_selenzyme.py "CC(=O)O>>CCO" --repo-root /path/to/selenzyme --host 511145 --targets 100
        """
    )
    
    parser.add_argument('smarts', help='SMARTS pattern for the reaction')
    parser.add_argument('--repo-root', help='Path to local SelenzymeRF repository')
    parser.add_argument('--api-url', default='http://selenzymerf.synbiochem.co.uk',
                       help='SelenzymeRF API URL')
    parser.add_argument('--targets', '-t', type=int, default=50,
                       help='Maximum number of results (default: 50)')
    parser.add_argument('--host', type=int, default=37919,
                       help='Host organism taxonomic ID (default: 37919 for R. opacus)')
    parser.add_argument('--output', '-o', help='Save results to CSV file')
    parser.add_argument('--no-msa', action='store_true',
                       help='Skip multiple sequence alignment')

    
    args = parser.parse_args()
    
    # Initialize client
    client = HybridSelenzymeRF(
        api_url=args.api_url,
        repo_root=args.repo_root
    )
    
    try:
        # Run search
        start_time = time.time()
        df = client.search_enzymes(
            smarts=args.smarts,
            targets=args.targets,
            host_taxid=args.host,
            nomsa=args.no_msa
        )
        elapsed = time.time() - start_time
        
        # Display results
        if len(df) > 0:
            print(f"\n🎉 Found {len(df)} enzyme suggestions in {elapsed:.1f}s")
            print(f"📊 Results source: {df['source'].iloc[0]}")
            
            # Show relevant columns
            display_cols = ['Score', 'Seq. ID', 'Description', 'Organism Source', 'Tax. distance', 'EC Number', 'Uniprot protein evidence', 'Consv. Score', 'Combined score']
            
            if display_cols:
                print("\nTop results:")
                print("=" * 80)
                print(df[display_cols].head(10).to_string(index=False))
            else:
                print(df.head(10).to_string(index=False))
            
            # Save if requested
            if args.output:
                df.to_csv(args.output, index=False)
                print(f"\n💾 Results saved to: {args.output}")
        else:
            print("❌ No enzyme suggestions found")
    
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()