import os
import pandas as pd
from ./classes/chem_utils import *
from ./classes/classes_functions import *
from ./classes/rp_classes import *

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
        if end_ik not in G.nodes:
            raise Exception("Warning: CAPE uses RetroPath only to connect the pollutant to Rhodococcus opacus PD630's central metabolites, but it does not connect the goal compound to them.")
        RP.django_call(name=start_name, inchi=start_inchi, steps=7)
	    RpOutput.tsv_generator(results_path="./cape_app/algorithms/pathways/rp_results/results.csv", outdir="./cape_app/algorithms/pathways/rp_results/")
	    Std.main(tsv_edges_path="./cape_app/algorithms/pathways/rp_results/rp_edges.tsv", edges_std_path="./cape_app/algorithms/pathways/rp_results/rp_std_edges.tsv", processes=64)
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

