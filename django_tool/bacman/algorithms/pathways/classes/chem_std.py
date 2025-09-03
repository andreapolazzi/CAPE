from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi
from rdkit import RDLogger

from multiprocessing import Pool
import functools

import pandas as pd
import numpy as np
import csv

# Classes for chemical standardization

class Std():

    @staticmethod
    def standardize(inchi_str, 
                standardize=True,
                neutralize=True, 
                strip_stereo=False, 
                canonicalize_tautomers=True):
        """
        Standardize molecule to KEGG-like representation.
        """
        RDLogger.DisableLog('rdApp.*')  # Disable RDKit warnings

        if inchi_str is None or pd.isna(inchi_str):
            return np.nan

        # Load molecule from InChI
        mol = Chem.MolFromInchi(inchi_str)
        if mol is None:
            return None
        
        #0. Standardize to avoid invalid groups
        if standardize:
            mol = rdMolStandardize.Normalize(mol)
        
        # 1. Neutralize charges if requested
        if neutralize:
            mol = rdMolStandardize.Uncharger().uncharge(mol)  
            #mol = uncharger.uncharge(mol)
        
        # 2. Canonicalize tautomers
        if canonicalize_tautomers:
            #mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
            try:
                mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
            except RuntimeError as e:
                # Likely caused by missing RingInfo
                try:
                    Chem.SanitizeMol(mol)
                    mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
                except Exception as inner_e:
                    #(f"Canonicalization failed for InChI: {inchi_str}\n{inner_e}")
                    return None
        
        # 3. Strip stereochemistry if requested
        if strip_stereo:
            Chem.RemoveStereochemistry(mol)
        
        # 4. Get standardized InChI and InChIKey
        std_inchi = inchi.MolToInchi(mol)
        std_inchikey = inchi.InchiToInchiKey(std_inchi)
        
        #return {"inchi": std_inchi, "inchikey": std_inchikey, "mol": mol}
        return std_inchikey
    
    @staticmethod
    def standardize_parallel(inchi_list, n_processes=32):
        """Parallelize standardization using multiprocessing"""
        standardize_func = functools.partial(Std.standardize)
        with Pool(processes=n_processes) as pool:
            results = pool.map(standardize_func, inchi_list)
        return results
    
    
class RpStd(Std):

    @staticmethod
    def _clean_rp_std(rpout_std_path):
        """
        Cleans Retropath results after standardization
        For each group of 'Substrate InChI' and 'Product InChI', retains all rows with the lowest 'Score'
        Assures no loss of reaction components
        Cleaning after standardization allows to retain better quality data than before standardization

        Parameters:
            rpout_std_path: path to standardized retropath output
        """

        df = pd.read_csv(rpout_std_path)

        df = df.dropna()
        
        df_sorted = df.sort_values(by='Score')  # sort by Score for readability

        min_scores = df_sorted.groupby(['std_inchikey_substrate', 'std_inchikey_product'])['Score'].transform('min')

        retained_df = df_sorted[min_scores == df_sorted['Score']]
        discarded_df = df_sorted[~df_sorted.index.isin(retained_df.index)]

        #print(f"Retained {len(retained_df)} rows, discarded {len(discarded_df)} rows")
        #print(retained_df, '\n', discarded_df)

        # For each unique 'Transformation ID' in retained_df, append matching rows from discarded_df
        # with the same 'Transformation ID'. This ensures all components of a reaction are retained
        transformation_ids = retained_df['Transformation ID'].unique()
        matching_discarded = discarded_df[discarded_df['Transformation ID'].isin(transformation_ids)]
        #print(f"Matching discarded rows: {len(matching_discarded)}")
        retained_df = pd.concat([retained_df, matching_discarded], ignore_index=True)

        #print(f"Taking back reaction components\n\
            #Retained {len(retained_df)} rows, discarded {len(discarded_df)} rows")
        #print(retained_df, '\n', discarded_df)

        # return retained_df, discarded_df
        return retained_df
    
    @staticmethod
    def rp_std(rp_res_path, rp_std_path):
        df = pd.read_csv(rp_res_path)

        #start_time = time.time()
        df['std_inchikey_substrate'] = Std.standardize_parallel(df['Substrate InChI'].tolist(), n_processes=32)
        df['std_inchikey_product'] = Std.standardize_parallel(df['Product InChI'].tolist(), n_processes=32)
        #end_time = time.time()
        df.to_csv(rp_std_path, index=False, quoting=csv.QUOTE_ALL)
        #print(f"Standardization completed in {end_time - start_time:.2f} seconds")

    @staticmethod
    def main(rp_res_path, rp_std_path):
                
        RpStd.rp_std(rp_res_path, rp_std_path)
        
        df = pd.read_csv(rp_std_path)
        #print(len(df), "rows processed")

        # Clean the results after standardization
        #ret_df, dis_df = clean_rpout(df)
        ret_df = RpStd._clean_rp_std(df)

        #return ret_df
        ret_df.to_csv(rp_std_path, index=False, quoting=csv.QUOTE_ALL)


class KeggStd(Std):
    
    @staticmethod
    def main(kegg_edges_path, kegg_std_path):
        """
        Standardizes KEGG edges
        Parameters:
            kegg_edges_path: path to KEGG edges file (tsv)
            kegg_std_path: path to save standardized KEGG edges (csv)
        """
    
        #df = pd.read_csv('/work/igem/standardization/roa_inchi_edges.tsv', sep="\t")
        df = pd.read_csv(kegg_edges_path, sep="\t")
        
        df.dropna(subset=['source_inchi', 'target_inchi'], inplace=True)
        # TODO: avoid having to dropna

        #start_time = time.time()
        df['std_inchikey_substrate'] = Std.standardize_parallel(df['source_inchi'].tolist(), n_processes=64)
        df['std_inchikey_product'] = Std.standardize_parallel(df['target_inchi'].tolist(), n_processes=64)
        #end_time = time.time()
        #df.to_csv('/work/igem/standardization/roa_std_edges.csv', index=False, quoting=csv.QUOTE_ALL)
        df.to_csv(kegg_std_path, index=False, quoting=csv.QUOTE_ALL)
        #print(df)
        #print(df.isna().sum())
        #print(f"Standardization completed in {end_time - start_time:.2f} seconds")