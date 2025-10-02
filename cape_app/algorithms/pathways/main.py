import os
import pandas as pd
from .chem_utils import *
from .classes_functions import *


def main(start_inchi: str, end_inchi: str, start_name:str='', end_name:str=''):
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

    # Check if start or end is missing in G
    if start_ik in G.nodes and end_ik in G.nodes:
        paths = k_shortest_paths(G, start_ik, end_ik, ik_to_name, K=MAX_PATHS)
        if paths: # Paths found, return the results
            return paths,info_dict

    # RetroPath fallback
    try:
        try_retropath(c_inchi=start_ik,c_name='')
        # Add new edges to G and update dictionaries with new info
        G, name_to_ik, ik_to_name, info_dict = build_graph_from_tsv(G,method="RetroPath")
        # Re-run path search
        paths = k_shortest_paths(G, start_ik, end_ik, ik_to_name, name_to_ik, K=MAX_PATHS)
        if not paths:
            return False # Still no paths found, even after successfully adding Retropath edges
        else: # Paths found, return the results
            return paths,info_dict
    except: # RetroPath failed
        print('retropath failed')
        return False
