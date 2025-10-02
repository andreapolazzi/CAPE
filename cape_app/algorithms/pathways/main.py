import os
from chem_utils import *
from classes_functions import *
import pandas as pd

def main(start_inchi: str, end_inchi: str):
    PATHWAY_FOLDER = 'Pathways'
    os.makedirs(PATHWAY_FOLDER, exist_ok=True)

    MAX_PATHS = 3

    if not os.listdir(PATHWAY_FOLDER):
        raise('Warning: the Pathways folder is missing or empty.'
            'Download pathways before attempting to find shortest paths.')

    G, name_to_ik, ik_to_name, info_dict = build_graph_from_tsv(method="KEGG")

    # Normalize inchis and standardize them to inchikeys
    start_ik = Std.standardize(start_inchi.strip())
    end_ik = Std.standardize(end_inchi.strip())

    # 1) Check if start or end is missing in G
    if start_ik in G.nodes and end_ik in G.nodes:
        paths = k_shortest_paths(G, start_ik, end_ik, ik_to_name, K=MAX_PATHS)
        if paths:
            return paths,info_dict # Paths found, return the results

    # RetroPath fallback
    try:
        try_retropath(c_inchi=start_ik,c_name='')
        # Add new edges to G and update dictionaries with new info
        G, name_to_ik, ik_to_name, info_dict = build_graph_from_tsv(G,method="RetroPath")
        # Re-run path search
        paths = k_shortest_paths(G, start_ik, end_ik, ik_to_name, name_to_ik, K=MAX_PATHS)
        if not paths:
            # Still no paths found, even after successfully adding Retropath edges
            return False
        else:
            # Paths found, return the results
            return paths,info_dict
    except:
        # RetroPath failed
        print('retropath failed')
        return False
    
print(main(start_inchi="InChI=1S/C8H6O4/c9-7(10)5-1-2-6(4-3-5)8(11)12/h1-4H,(H,9,10)(H,11,12)",end_inchi="InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)")[0])
# terephthalate InChI=1S/C8H6O4/c9-7(10)5-1-2-6(4-3-5)8(11)12/h1-4H,(H,9,10)(H,11,12)
# to pyruvate