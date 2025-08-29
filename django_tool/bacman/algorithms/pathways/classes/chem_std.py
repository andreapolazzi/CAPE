from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi
from rdkit import RDLogger

from multiprocessing import Pool
import functools

import pandas as pd

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
    
    @staticmethod
    def clean_rpout(rpout_std):
        """
        Cleans Retropath results after standardization
        For each group of 'Substrate InChI' and 'Product InChI', retains all rows with the lowest 'Score'
        Assures no loss of reaction components
        Cleaning after standardization allows to retain better quality data than before standardization

        Parameters:
            rpout_std: path to standardized retropath output
        """

        df = pd.read_csv(rpout_std)

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

        return retained_df, discarded_df