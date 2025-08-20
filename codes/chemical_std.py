from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi
from rdkit import RDLogger

from multiprocessing import Pool
import functools

import pandas as pd
import time
import csv
import matplotlib.pyplot as plt

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
        #tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
        #mol = tautomer_enumerator.Canonicalize(mol)
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    
    # 3. Strip stereochemistry if requested
    if strip_stereo:
        Chem.RemoveStereochemistry(mol)
    
    # 4. Get standardized InChI and InChIKey
    std_inchi = inchi.MolToInchi(mol)
    std_inchikey = inchi.InchiToInchiKey(std_inchi)
    
    #return {"inchi": std_inchi, "inchikey": std_inchikey, "mol": mol}
    return std_inchikey

def standardize_parallel(inchi_list, n_processes=2):
    """Parallelize standardization using multiprocessing"""
    standardize_func = functools.partial(standardize)
    with Pool(processes=n_processes) as pool:
        results = pool.map(standardize_func, inchi_list)
    return results

def clean_rpout(df):
    """
    Cleans Retropath results after standardization
    For each group of 'Substrate InChI' and 'Product InChI', retains all rows with the lowest 'Score'
    Assures no loss of reaction components
    Cleaning after standardization allows to retain better quality data than before standardization
    """

    df = df.dropna()
    
    df_sorted = df.sort_values(by='Score')  # sort by Score for readability

    min_scores = df_sorted.groupby(['std_inchikey_substrate', 'std_inchikey_product'])['Score'].transform('min')

    retained_df = df_sorted[min_scores == df_sorted['Score']]
    discarded_df = df_sorted[~df_sorted.index.isin(retained_df.index)]

    print(f"Retained {len(retained_df)} rows, discarded {len(discarded_df)} rows")
    print(retained_df, '\n', discarded_df)

    # For each unique 'Transformation ID' in retained_df, append matching rows from discarded_df
    # with the same 'Transformation ID'. This ensures all components of a reaction are retained
    transformation_ids = retained_df['Transformation ID'].unique()
    matching_discarded = discarded_df[discarded_df['Transformation ID'].isin(transformation_ids)]
    print(f"Matching discarded rows: {len(matching_discarded)}")
    retained_df = pd.concat([retained_df, matching_discarded], ignore_index=True)

    print(f"Taking back reaction components\n\
          Retained {len(retained_df)} rows, discarded {len(discarded_df)} rows")
    print(retained_df, '\n', discarded_df)

    return retained_df, discarded_df

def boxplots(retained_df, discarded_df, path_fname):
        """
        Plots two boxplots comparing the 'Score' distributions of retained and discarded DataFrames.
        """
        scores_data = [retained_df['Score'], discarded_df['Score']]
        labels = ['Retained', 'Discarded']

        plt.figure(figsize=(8, 5))
        plt.boxplot(scores_data, tick_labels=labels)
        plt.title('Score Distribution: Retained vs Discarded')
        plt.ylabel('Score')
        plt.grid(True)
        plt.savefig(path_fname)

def main_bp3():

    RDLogger.DisableLog('rdApp.*')  # Disable RDKit warnings
    """
    df = pd.read_csv('/work/igem/iGEM/results/BP3_6_8_100_7steps/results.csv')

    start_time = time.time()
    df['std_inchikey_substrate'] = standardize_parallel(df['Substrate InChI'].tolist(), n_processes=128)
    df['std_inchikey_product'] = standardize_parallel(df['Product InChI'].tolist(), n_processes=128)
    end_time = time.time()
    df.to_csv('/work/igem/iGEM/results/BP3_6_8_100_7steps/results_std.csv', index=False, quoting=csv.QUOTE_ALL)
    print(f"Standardization completed in {end_time - start_time:.2f} seconds")
    """
    df = pd.read_csv('/work/igem/iGEM/results/BP3_6_8_100_7steps/results_std.csv')
    print(len(df), "rows processed")

    # Clean the results after standardization
    ret_df, dis_df = clean_rpout(df)
    boxplots(ret_df, dis_df, '/work/igem/iGEM/results/BP3_6_8_100_7steps/ret_disc_std.png')
    # after standardization, managed to clean 7000 rows
    # (instead of 6500 doing it before standardization)
    # but regained more, so both methods produce about 21500 retained rows
    # still good because when there will be more duplicates, might not regain so many
    # BUT if dropna() is used, there is almost no difference between two methods

if __name__ == '__main__':

    main_bp3()
    # TODO: Caching

    """
    # Example: L-lactic acid InChI (charged vs neutral)
    lactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)"
    lactate_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1"
    Dlactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m1/s1"
    Llactic_inchi = "InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1"
    pyruvate_inchi = "InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)/p-1"

    # standardize it in KEGG-like style
    standardized_lactic = standardize(lactic_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)
    standardized_lactate = standardize(lactate_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)
    standardized_Dlactic = standardize(Dlactic_inchi, 
                                    neutralize=True, 
                                    strip_stereo=False, 
                                    canonicalize_tautomers=False)

    print("standardized lactic InChI:", standardized_lactic["inchi"])
    print("standardized lactic InChIKey:", standardized_lactic["inchikey"])
    print("standardized lactate InChI:", standardized_lactate["inchi"])
    print("standardized lactate InChIKey:", standardized_lactate["inchikey"])
    print("standardized pyruvate InChI:", standardized_Dlactic["inchi"])
    print("standardized pyruvate InChIKey:", standardized_Dlactic["inchikey"])

    m = Chem.MolFromInchi(Llactic_inchi)
    m2 = rdMolStandardize.Normalize(m)
    m2_inchi = inchi.MolToInchi(m2)
    m2_inchikey = inchi.InchiToInchiKey(m2_inchi)
    print(m2_inchi, m2_inchikey)
    """
