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
    def standardize(inchi_str: str, 
                standardize=True,
                neutralize=True, 
                strip_stereo=False, 
                canonicalize_tautomers=True) -> str:
        """
        Standardizes molecule to KEGG-like representation.
        Args:
            inchi_str (str): InChI string of the molecule to standardize.
            standardize (bool): Whether to apply normalization.
            neutralize (bool): Whether to neutralize charges.
            strip_stereo (bool): Whether to remove stereochemistry.
            canonicalize_tautomers (bool): Whether to canonicalize tautomers.
        Returns:
            str: Standardized InChIKey or None if standardization fails.
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
    
    @staticmethod
    def standardize_parallel(inchi_list: list, n_processes=2) -> list:
        """
        Parallelize standardization using multiprocessing
        Args:
            inchi_list (list): List of InChI strings to standardize.
            n_processes (int): Number of parallel processes to use.
        Returns:
            results: List of standardized InChIKeys or None for failed standardizations.
        """
        standardize_func = functools.partial(Std.standardize)
        with Pool(processes=n_processes) as pool:
            results = pool.map(standardize_func, inchi_list)
        return results
    
    @staticmethod
    def std(tsv_edges_path: str, n_processes=2) -> pd.DataFrame:
        """
        Standardizes compounds in the edges file.
        Args:
            tsv_edges_path (str): Path to the TSV file containing edges 
                with 'source_inchi' and 'target_inchi' columns.
            n_processes (int): Number of parallel processes to use.
        Returns:
            pd.DataFrame: DataFrame with added 
                'source_inchi_key' and 'target_inchi_key' columns.
        """
        df = pd.read_csv(tsv_edges_path, sep="\t")
        df['source_inchi_key'] = Std.standardize_parallel(df['source_inchi'].tolist(), 
                                                          n_processes=n_processes)
        df['target_inchi_key'] = Std.standardize_parallel(df['target_inchi'].tolist(), 
                                                          n_processes=n_processes)
        return df
    
    @staticmethod
    def main(tsv_edges_path: str, edges_std_path:str, processes=2, silent=True):
        """
        Args:
            tsv_edges_path: path to tsv edges file to standardize
                the path must end with a .tsv -> example: /path/to/rp_output.tsv
            edges_std_path: path to save standardized edges output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_std_output.tsv
        """
        if silent:
            RDLogger.DisableLog('rdApp.*')  # Disable RDKit warnings
        df = Std.std(tsv_edges_path, n_processes=processes)  # standardizes nodes
        df.dropna(inplace=True)  # drop rows where standardization failed
        df.to_csv(edges_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)
    
    
class RpStd(Std):

    @staticmethod
    def _clean_rp_std(rpout_std_path: str) -> pd.DataFrame:
        """
        Cleans Retropath results after standardization
        For each group of 'Substrate InChI' and 'Product InChI', retains all rows with the lowest 'Score'
        Assures no loss of reaction components
        Cleaning after standardization allows to retain better quality data than before standardization

        Args:
            rpout_std_path: path to standardized retropath output
        Returns:
            retained_df: dataframe with cleaned retropath results
        """
        df = pd.read_csv(rpout_std_path)
        df = df.dropna()
        df_sorted = df.sort_values(by='Score')  # sort by Score for readability
        min_scores = df_sorted.groupby(['std_inchikey_substrate', 'std_inchikey_product'])['Score'].transform('min')

        retained_df = df_sorted[min_scores == df_sorted['Score']]
        discarded_df = df_sorted[~df_sorted.index.isin(retained_df.index)]

        # For each unique 'Transformation ID' in retained_df, append matching rows from discarded_df
        # with the same 'Transformation ID'. This ensures all components of a reaction are retained
        transformation_ids = retained_df['Transformation ID'].unique()
        matching_discarded = discarded_df[discarded_df['Transformation ID'].isin(transformation_ids)]
        retained_df = pd.concat([retained_df, matching_discarded], ignore_index=True)

        return retained_df

    @staticmethod
    def main(rp_res_path: str, rp_std_path: str, keep_best=False):
        """
        Main for standardizing specifically retropath results.
        Args:
            rp_res_path: path to retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_output.tsv
            rp_std_path: path to save standardized retropath output (tsv)
                the path must end with a .tsv -> example: /path/to/rp_std_output.tsv
            keep_best: whether to clean the results after standardization
        Returns:
            None
        """
        df = Std.std(rp_res_path)
        
        if keep_best:
            work_df = pd.read_csv(rp_std_path)
            ret_df = RpStd._clean_rp_std(work_df)
            df = ret_df.copy(deep=True)

        df.dropna(inplace=True)
        df.to_csv(rp_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)


class KeggStd(Std):
    
    @staticmethod
    def main(kegg_edges_path: str, kegg_std_path: str):
        """
        Standardizes KEGG edges.
        Warning: this function is not up to date with the rest of the codebase
        Args:
            kegg_edges_path: path to KEGG edges file (tsv)
            kegg_std_path: path to save standardized KEGG edges (csv)
        """
        df = pd.read_csv(kegg_edges_path, sep="\t")
        df.dropna(subset=['source_inchi', 'target_inchi'], inplace=True)
        # TODO: avoid having to dropna

        df['std_inchikey_substrate'] = Std.standardize_parallel(df['source_inchi'].tolist(), n_processes=64)
        df['std_inchikey_product'] = Std.standardize_parallel(df['target_inchi'].tolist(), n_processes=64)
        df.to_csv(kegg_std_path, sep='\t', index=False, quoting=csv.QUOTE_ALL)


class Penalty():

    @staticmethod
    def transf_distance(substrate_id: str, product_id: str, format="inchi") -> float:
        """
        Uses Tanimoto similarity of Morgan fingerprints to estimate chemical distance
        Args:
            substrate_id: SMILES or InChI of the substrate
            product_id: SMILES or InChI of the product
            format: "smiles" or "inchi" (default)
        Returns:
            1 - similarity: distance (higher = more chemical change)
        """

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
    def dist_penalty(distance: float) -> float:
        """
        Converts chemical distance to a penalty score to be used as edge weight.
        Highly penalizes large chemical changes, e.g. from a metabolite to a cofactor.
        Args:
            distance: distance from transf_distance() (0-1)
        Returns:
            penalty: 3 for highly probable cofactors or currency metabolites
        """
        # if distance is not None:
        #     if distance <= 0.70:
        #         # Main pathway
        #         penalty = distance * 2.0
        #     elif distance < 0.80:
        #         # Still main pathway - some doubts
        #         penalty = 0.60 * 2.0 + (distance - 0.60) * 5.0
        #     else:
        #         # Cofactors or byproducts
        #         penalty = 0.60 * 2.0 + 0.20 * 5.0 + (distance - 0.80) * 12.0  # max add 4.6

        #     return penalty
        # else:
        #     return None
        if distance is not None:
            if distance < 0.99:
                penalty = 0
            else:
                # Cofactors or byproducts
                penalty = (distance - 0.99) * 300
            return penalty
        else:
            return None
        
    @staticmethod
    def get_mw(identifier: str, format="inchi") -> float:
        """
        Calculates molecular weight from SMILES or InChI
        Args:
            identifier: SMILES or InChI of the molecule
            format: "smiles" or "inchi" (default)
        Returns:
            Molecular weight or None if calculation fails
        """
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
    def mw_change(substrate_mw: float, product_mw: float) -> float:
        """
        Molecular Weight-based penalty
        Args:
            substrate_mw: molecular weight of the substrate
            product_mw: molecular weight of the product
        Returns:
            penalty: molecular weight change penalty (0-3)
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
    def base_penalty(color: str) -> float:
        """
        Adds penalty to favor:
            green kegg edges (annotated in PD630), over
            white kegg edges (not annotated in PD630), over
            RetroPath edges (predicted reactions)
        Args:
            color: "green", "white", "RetroPath"
        Returns:
            penalty: 0 for green, 2 for white, 4 for RetroPath
        """
        if color == "green":
            return 0.0
        elif color == "white":
            return 2.0
        elif color == "RetroPath":
            return 4.0
        else:
            return 10  # default for unknown colors
        
    @staticmethod
    def transf_penalty(base_penalty: float, distance_penalty: float, mw_penalty: float, rp_score: float) -> float:
        """
        Penalty Assembler
        
        Args:
            base_penalty: discrimnant for kegg/retropath edges
            distance_penalty: adjusted Tanimoto distance -> 0-4.6 (0=identical, 4.6=completely different)
            mw_penalty: molecular weight change -> 0-3 (0=no change, 3=huge mw change)
            rp_score: Retropath score -> 0-4 (0=best, 4=worst)
        
        Returns:
            Combined penalty score (weight to be used as edge weight)
        """
        if distance_penalty is not None and mw_penalty is not None:
            penalty = base_penalty + distance_penalty + mw_penalty + rp_score
            return penalty
        else:
            return None
        
    @staticmethod
    def main(std_path: str, weighted_path: str, write=True) -> pd.DataFrame:
        """
        Generates the weighted dataframe needed for the weighted graph

        Args:
            std_path: path to standardized edges tsv file
            weighted_path: path to save the weighted edges tsv file
        Returns:
            weighted_df: dataframe with weights
        """
        
        # std_df = pd.read_csv(std_path, sep="\t")
        # invalid_df = pd.read_csv("./cape_app/algorithms/pathways/data/invalid_nodes.tsv", sep="\t")
        # invalid_iks = set(invalid_df['std_inchikey'].tolist())

        # Create mask for invalid rows
        # is_invalid = std_df['source_inchi_key'].isin(invalid_iks) | std_df['target_inchi_key'].isin(invalid_iks)
        # valid_mask = ~is_invalid
        
        # # Only compute distances and penalties for valid rows
        # std_df.loc[valid_mask, "distance"] = std_df[valid_mask].apply(
        #     lambda row: Penalty.transf_distance(row["source_inchi"], row["target_inchi"]), axis=1
        # )
        # std_df.loc[valid_mask, "dist_penalty"] = std_df[valid_mask].apply(
        #     lambda row: Penalty.dist_penalty(row["distance"]), axis=1
        # )
        # std_df.loc[valid_mask, "mw_change"] = std_df[valid_mask].apply(
        #     lambda row: Penalty.mw_change(Penalty.get_mw(row["source_inchi"]), Penalty.get_mw(row["target_inchi"])), axis=1
        # )

        # weighted_df = std_df.copy(deep=True)
        # weighted_df.loc[is_invalid, "weight"] = 1000000  # Assign high weight to invalid rows
        # # Calculate normal weight for valid rows
        # weighted_df.loc[valid_mask, "weight"] = weighted_df[valid_mask].apply(
        #     lambda row: Penalty.transf_penalty(row["RPscore"], row["dist_penalty"], row["mw_change"]), axis=1
        # )
        # weighted_df["weight"] = weighted_df["weight"].fillna(20)

        # weighted_df.to_csv(weighted_path, sep="\t", index=False, quoting=csv.QUOTE_ALL)
        # return weighted_df

        std_df = pd.read_csv(std_path, sep="\t")
        invalid_df = pd.read_csv("/work/igem/iGEM/cape_app/algorithms/pathways/data/invalid_nodes.tsv", sep="\t")
        invalid_iks = set(invalid_df['std_inchikey'].tolist())

        # Create mask for invalid rows
        is_invalid = std_df['source_inchi_key'].isin(invalid_iks) | std_df['target_inchi_key'].isin(invalid_iks)
        valid_mask = ~is_invalid
        
        # Only compute distances and penalties for valid rows
        std_df.loc[valid_mask, "distance"] = std_df[valid_mask].apply(
            lambda row: Penalty.transf_distance(row["source_inchi"], row["target_inchi"]), axis=1
        )
        std_df.loc[valid_mask, "dist_penalty"] = std_df[valid_mask].apply(
            lambda row: Penalty.dist_penalty(row["distance"]), axis=1
        )
        std_df.loc[valid_mask, "mw_change"] = std_df[valid_mask].apply(
            lambda row: Penalty.mw_change(Penalty.get_mw(row["source_inchi"]), Penalty.get_mw(row["target_inchi"])), axis=1
        )
        std_df.loc[valid_mask, "base_penalty"] = std_df[valid_mask].apply(
            lambda row: Penalty.base_penalty(row["color"]), axis=1
        )
        
        weighted_df = std_df.copy(deep=True)
        weighted_df.loc[is_invalid, "weight"] = 1000000  # Assign high weight to invalid rows
        # Calculate normal weight for valid rows
        weighted_df.loc[valid_mask, "weight"] = weighted_df[valid_mask].apply(
            lambda row: Penalty.transf_penalty(row["base_penalty"], row["dist_penalty"], 
                                               row["mw_change"], row.get("Score", 0)), axis=1
        )
        weighted_df["weight"] = weighted_df["weight"].fillna(50)

        if write:
            weighted_df.to_csv(weighted_path, sep="\t", index=False, quoting=csv.QUOTE_ALL)
        return weighted_df