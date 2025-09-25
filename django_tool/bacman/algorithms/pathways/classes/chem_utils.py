from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors

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
        
        if inchi_str is None or pd.isna(inchi_str):
            return None

        try:
            mol = Chem.MolFromInchi(inchi_str)
            if mol is None:
                raise ValueError("MolFromInchi failed")

            if standardize:
                normalizer = rdMolStandardize.Normalizer()
                mol = normalizer.normalize(mol)

            if neutralize:
                uncharger = rdMolStandardize.Uncharger()
                mol = uncharger.uncharge(mol)

            if canonicalize_tautomers:
                enumerator = rdMolStandardize.TautomerEnumerator()
                mol = enumerator.Canonicalize(mol)

            if strip_stereo:
                Chem.RemoveStereochemistry(mol)

            Chem.SanitizeMol(mol)

            std_inchi = inchi.MolToInchi(mol)
            std_inchikey = inchi.InchiToInchiKey(std_inchi)

            return std_inchikey

        except Exception as e:
            return None
        
        """ old version
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
        """
    
    @staticmethod
    def standardize_parallel(inchi_list, n_processes=2):
        """Parallelize standardization using multiprocessing"""
        standardize_func = functools.partial(Std.standardize)
        with Pool(processes=n_processes) as pool:
            results = pool.map(standardize_func, inchi_list)
        return results
    
    @staticmethod
    def std(tsv_edges_path):
        df = pd.read_csv(tsv_edges_path, sep="\t")

        #start_time = time.time()
        df['source_inchi_key'] = Std.standardize_parallel(df['source_inchi'].tolist(), n_processes=2)
        df['target_inchi_key'] = Std.standardize_parallel(df['target_inchi'].tolist(), n_processes=2)
        #end_time = time.time()
        #print(f"Standardization completed in {end_time - start_time:.2f} seconds")
        return df
    
    @staticmethod
    def main(tsv_edges_path, edges_std_path):
        """
        Args:
            rp_res_path: path to retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_output.tsv
            rp_std_path: path to save standardized retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_std_output.tsv
        """
        df = Std.std(tsv_edges_path)

        df.dropna(inplace=True)

        df.to_csv(edges_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)
    
    
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
    def main(rp_res_path, rp_std_path, keep_best=False):
        """
        Args:
            rp_res_path: path to retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_output.tsv
            rp_std_path: path to save standardized retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_std_output.tsv
        """
        df = Std.std(rp_res_path)
        
        if keep_best:
            work_df = pd.read_csv(rp_std_path)
            #print(len(df), "rows processed")

            # Clean the results after standardization
            #ret_df, dis_df = clean_rpout(df)
            ret_df = RpStd._clean_rp_std(work_df)
            df = ret_df.copy(deep=True)

        #return ret_df
        df.dropna(inplace=True)
        df.to_csv(rp_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)


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
        df.to_csv(kegg_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)
        #print(df)
        #print(df.isna().sum())
        #print(f"Standardization completed in {end_time - start_time:.2f} seconds")



class Penalty():

    @staticmethod
    def transf_distance(substrate_id, product_id, format="inchi"):

        if format == "smiles":
            mol1 = Chem.MolFromSmiles(substrate_id)
            mol2 = Chem.MolFromSmiles(product_id)
        else:  # default InChI
            mol1 = Chem.MolFromInchi(substrate_id)
            mol2 = Chem.MolFromInchi(product_id)
        
        if mol1 is None or mol2 is None:
            return None  # Invalid SMILES/InChI

        try:
            # Generate Morgan fingerprints
            fp1 = AllChem.GetMorganFingerprint(mol1, 2)  # radius=2
            fp2 = AllChem.GetMorganFingerprint(mol2, 2)
            
            # Calculate similarity
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            
            # Return transformation significance (1 - similarity)
            return 1 - similarity  # Higher = more chemical change
        except:
            return None  # Error in processing
    
    @staticmethod
    def dist_penalty(distance):
        if distance is not None:
            if distance <= 0.70:
                # Main pathway
                penalty = distance * 2.0
                #penalty += 1.4
            elif distance < 0.80:
                # Still main pathway - some doubts
                penalty = 0.60 * 2.0 + (distance - 0.60) * 5.0
            else:
                # Cofactors or byproducts
                penalty = 0.60 * 2.0 + 0.20 * 5.0 + (distance - 0.80) * 12.0  # max add 4.6
        
            return penalty
        else:
            return None
        
    @staticmethod
    def get_mw(identifier, format="inchi"):
        try:
            if format == "smiles":
                mol = Chem.MolFromSmiles(identifier)
            else:  # default InChI
                mol = Chem.MolFromInchi(identifier)
            if mol is not None:
                return Descriptors.MolWt(mol)
            else:
                return None
        except:
            return None
    
    @staticmethod
    def mw_change(substrate_mw, product_mw):
        """
        Molecular Weight based penalty
        """
        if substrate_mw is not None and product_mw is not None and substrate_mw > 0:
            mw_change = substrate_mw - product_mw  # Positive = loss, Negative = gain
            percent_change = abs(mw_change) / substrate_mw * 100
            
            if mw_change >= 0:  # MW LOSS
                # Exponential curve: gentle rise then steep
                penalty = 3.0 * (1 - np.exp(-percent_change / 40))
                return min(penalty, 3.0)
            else:  # MW GAIN - more aggressive penalty
                # Steeper exponential for gains
                penalty = 3.0 * (1 - np.exp(-percent_change / 25))
                return min(penalty, 3.0)
        else:
            return None
        
    @staticmethod
    def transf_penalty(rp_score, distance_penalty, mw_penalty):
        """
        Penalty Assembler
        
        Args:
            rp_score: RetroPath reaction score -> 0-4 (0=best, 4=worst)
                Use 0 for KEGG reactions
            distance_penalty: adjusted Tanimoto distance -> 0-4.6 (0=identical, 4.6=completely different)
            mw_penalty: molecular weight change -> 0-3 (0=no change, 3=huge mw change)
            kegg_penalty: penalty for non-KEGG reactions -> default 2
                Use 0 for KEGG reactions
        
        Returns:
            Combined penalty score (weight to be used as edge weight)
        """
        if distance_penalty is not None and mw_penalty is not None:
            penalty = rp_score + distance_penalty + mw_penalty
            return penalty
        else:
            return None
        
    @staticmethod
    def main(std_df):
        """
        Generates the weighted dataframe needed for the weighted graph

        Args:
            std_df: dataframe from KEGG or RetroPath
        Returns:
            weighted_df: dataframe with weights
        """
        std_df["distance"] = std_df.apply(lambda row: Penalty.transf_distance(row["source_inchi"], row["target_inchi"]), axis=1)
        std_df["dist_penalty"] = std_df.apply(lambda row: Penalty.dist_penalty(row["distance"]), axis=1)
        std_df["mw_change"] = std_df.apply(lambda row: Penalty.mw_change(Penalty.get_mw(row["source_inchi"]), Penalty.get_mw(row["target_inchi"])), axis=1)

        weighted_df = std_df.copy(deep=True)
        weighted_df["weight"] = weighted_df.apply(lambda row: Penalty.transf_penalty(row["RPscore"], row["dist_penalty"], row["mw_change"]), axis=1)
        
        return weighted_df