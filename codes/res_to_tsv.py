import pandas as pd
from pubchempy import get_compounds
import requests

def load_res(file_path):
    """
    Load the results from a CSV file into a pandas DataFrame.
    """
    return pd.read_csv(file_path)

def get_name(inchi):
        """
        Get the name of a compound from its InChI ID.
        """
        if inchi == 'InChI=1S':
            return 'H+'
        compounds = get_compounds(inchi, 'inchi')
        if compounds:
            if compounds[0].synonyms:
                return compounds[0].synonyms[0]
            elif compounds[0].iupac_name:
                return compounds[0].iupac_name
            else:
                return 'Unknown'
        return 'Unknown'

def add_names(df):
    """
    Add 'Substrate Name' and 'Product Name' columns to the DataFrame based on InChI IDs
    Queries pubchem to get names.
    """

    df['Substrate Name'] = df['Substrate InChI'].apply(get_name)
    df['Product Name'] = df['Product InChI'].apply(get_name)
    
    return df

def nodes_to_tsv(df, outdir):
    """
    Save the unique substrate and product names to a TSV file.
    """
    
    unique_substrates = pd.Series(df['Substrate InChI'].unique())
    unique_products = pd.Series(df['Product InChI'].unique())

    all_compounds = pd.concat([unique_substrates, unique_products]).unique()
    all_compounds = pd.DataFrame({'inchi': all_compounds})
    all_compounds['name'] = all_compounds['inchi'].apply(get_name)
    all_compounds = all_compounds[['name', 'inchi']]

    all_compounds.to_csv(outdir+'rp_nodes.tsv', sep='\t', index=False, header=True)


def ec_to_name(ec: str):
    """
    Convert EC number (format "x.x.x.x", not "EC:x.x.x.x") to a human-readable name.
    """
    try:
        response = requests.get(f'https://rest.kegg.jp/get/ec:{ec}')
        if response.status_code == 200:
            lines = response.text.split('\n')
            for line in lines:
                if line.startswith('NAME'):
                    return line.split('NAME')[1].strip().strip(';')
        return 'Unknown'
    except:
        return 'Unknown'


def edges_to_tsv(df, outdir):
    """
    Save the edges (reactions) to a TSV file.
    """
    df['EC number'] = df['EC number'].apply(lambda x: x.strip('[]'))
    #print(df['EC number'])
    edges_df = df[['EC number', 'Substrate InChI', 'Product InChI']]
    #print(edges_df)
    edges_df['enzyme_name'] = edges_df['EC number'].apply(lambda x: x if x == 'NOEC' else ec_to_name(x))
    edges_df['color'] = 'yellow'
    edges_df = edges_df.rename(columns={'EC number': 'EC', 'Substrate InChI': 'source', 'Product InChI': 'target'})
    edges_df = edges_df[['enzyme_name', 'EC', 'source', 'target', 'color']]
    edges_df.to_csv(outdir+'rp_edges.tsv', sep='\t', index=False, header=True)
    #print(edges_df)
    


if __name__ == "__main__":
    
    #file_path = "/work/igem/iGEM/data_files/res_cases/ota_MET_scope.csv"
    file_path = "./data_files/res_cases/ota_MET_scope.csv"
    outdir = "./results/"
    
    df = load_res(file_path)
    #print("Initial df\n", df)
    
    nodes_to_tsv(df, outdir)
    print("Nodes saved to TSV file.")

    #print(ec_to_name("1.18.1.3"))
    edges_to_tsv(df, outdir)
    print("Edges saved to TSV file.")
    
    


    