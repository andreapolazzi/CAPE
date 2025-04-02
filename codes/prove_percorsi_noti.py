from prova_launch import launch_retropath
import requests
from rdkit import Chem

def get_inchi_from_kegg(kegg_id):
    """
    Retrieve InChI for a molecule given its KEGG ID
    Args: kegg_id (str): KEGG compound ID (e.g., 'C00031')   
    Returns: str: InChI string or None if unsuccessful
    """
    # Get the MOL file from KEGG
    url = f"https://rest.kegg.jp/get/{kegg_id}/mol"
    response = requests.get(url)
    
    if response.status_code != 200 or not response.text.strip():
        print(f"Failed to retrieve MOL file for {kegg_id}")
        return None
    
    # Convert MOL to InChI using RDKit
    try:
        mol = Chem.MolFromMolBlock(response.text)
        if mol is None:
            print(f"Failed to parse MOL file for {kegg_id}")
            return None
        
        inchi = Chem.MolToInchi(mol)
        return inchi
    
    except Exception as e:
        print(f"Error converting to InChI: {e}")
        return None

def glu_f6p():
    kegg_ids = [
                "C00668",  # glucose 6-phosphate
                "C00085",  # fructose 6-phosphate
                "C00354",  # fructose 1,6-bisphosphate
                "C00118",  # glyceraldehyde 3-phosphate
                ]
    
    inchis = {kegg_id:get_inchi_from_kegg(kegg_id) for kegg_id in kegg_ids}
    if inchis:
        print(f"{inchis}")

    launch_retropath(rules_file="iGEM/data_files/noti/rules1.csv",
                     source_file="iGEM/data_files/noti/source1.csv",
                     sink_file="iGEM/data_files/noti/sink1.csv")


def get_retrorules_from_inchi(inchi, diameter=12):
    """
    Query RetroRules API to get reaction rules based on the provided InChI
    Args:
        inchi (str): InChI string of the compound
        diameter (int): The reaction rule diameter (default: 12)
    Returns:
        list: List of reaction rules or None if unsuccessful
    """
    try:
        inchi_url = "https://retrorules.org/api/v0.7/inchi/?"
        
        input_query = "input="
        minDiameter = "minDiameter="

        input_query += inchi
        print(input_query)
        minDiameter += str(diameter)
        print(minDiameter)

        url = inchi_url + input_query + "&" + minDiameter
        print(url)
        
        response = requests.get(url)
        
        if response.status_code != 200:
            print(f"API request failed with status code {response.status_code}")
            return None
        
        # Parse the JSON response
        data = response.json()
        return data.get("rules", [])
    
    except Exception as e:
        print(f"Error querying RetroRules API: {e}")
        return None

def glu_g3p():

    kegg_ids = [
                "C00668",  # glucose 6-phosphate
                "C00085",  # fructose 6-phosphate
                "C00354",  # fructose 1,6-bisphosphate
                "C00118",  # glyceraldehyde 3-phosphate
                ]
    inchis = {kegg_id:get_inchi_from_kegg(kegg_id) for kegg_id in kegg_ids}

    for kId, inchi in inchis.items():
        print(kId, inchi)
    #get_retrorules_from_inchi(inchis["C00668"])

    # insert code here

    launch_retropath(#rules_file="iGEM/data_files/noti/retrorules_rr02_rp2_flat_forward.csv",
                     rules_file="iGEM/data_files/noti/rules2.csv",
                     source_file="iGEM/data_files/noti/source2.csv",
                     sink_file="iGEM/data_files/noti/sink2.csv")


# Example usage
if __name__ == "__main__":

    #glu_f6p()
    glu_g3p()