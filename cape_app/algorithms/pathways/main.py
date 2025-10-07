import os
import pandas as pd
from .classes.chem_utils import *
from .classes.classes_functions import *
from .classes.rp_classes import *

def main(start_inchi: str, end_inchi: str, start_name:str='', end_name:str='', max_paths: int = 3,
         outdir="./cape_app/pipeline_results/rp_results/",
         max_steps=5, max_structures=100, min_diameter=6, max_diameter=16, timeout=60,
         sink_file="./cape_app/algorithms/pathways/data/gsm_gapfilled_sink.csv",
         rules_file="./cape_app/algorithms/pathways/data/retrorules_rr02_rp2_flat_forward.csv"):
    
    script_dir = os.path.dirname(os.path.abspath(__file__))

    PATHWAY_FOLDER = os.path.join(script_dir, 'data', 'Pathways')
    os.makedirs(PATHWAY_FOLDER, exist_ok=True)

    if max_paths <= 0:
        max_paths = 3

    if not os.listdir(PATHWAY_FOLDER):
        raise('Warning: the Pathways folder is missing or empty.'
            'Download pathways before attempting to find shortest paths.')

    G = build_graph_from_tsv(method="KEGG")

    # Normalize inchis and standardize them to inchikeys
    start_ik = Std.standardize(start_inchi.strip())
    end_ik = Std.standardize(end_inchi.strip())

    # Check if start or end is missing in G
    if start_ik in G.nodes and end_ik in G.nodes:
        print("trying kegg")
        paths = k_shortest_paths(G, start_ik, end_ik, K=max_paths)
        if paths: # Paths found, return the results
            return paths

    # RetroPath fallback
    try:
        if end_ik not in G.nodes:
            print(("Warning: CAPE uses RetroPath only to connect the pollutant " \
            "to Rhodococcus opacus PD630's central metabolites, but it " \
            "does not connect the goal compound to them."))
            raise Exception("Warning: CAPE uses RetroPath only to connect the pollutant " \
            "to Rhodococcus opacus PD630's central metabolites, but it " \
            "does not connect the goal compound to them.")
        
        print('trying to call retropath')
        # TODO: adjust with input outdir
        os.makedirs(outdir, exist_ok=True)
        RP.django_call(name=start_name, inchi=start_inchi, steps=max_steps,
                       sink=sink_file, rules=rules_file,
                       structures=max_structures, dmin=min_diameter, dmax=max_diameter,
                       outdir=outdir, timeout=timeout)
        print("generating tsv from retropath results")
        RpOutput.tsv_generator(results_path="./cape_app/pipeline_results/rp_results/results.csv", 
                               outdir="./cape_app/pipeline_results/rp_results/")
        print("standardizing retropath edges")
        Std.main(tsv_edges_path="./cape_app/pipeline_results/rp_results/rp_edges.tsv", 
                 edges_std_path="./cape_app/pipeline_results/rp_results/rp_std_edges.tsv", 
                 processes=64)
        print("adding weights to retropath edges")
        Penalty.main(std_path="./cape_app/pipeline_results/rp_results/rp_std_edges.tsv", 
                     weighted_path="./cape_app/pipeline_results/rp_results/rp_weighted_edges.tsv")
        
        # Add new edges to G and update dictionaries with new info
        print("trying to add retropath edges to graph")
        G = build_graph_from_tsv(G, method="RetroPath",
                                 rp_outdir="./cape_app/pipeline_results/rp_results/")
        # Re-run path search
        print("re-running path search")
        paths = k_shortest_paths(G, start_ik, end_ik, K=max_paths)
        if not paths:
            return False # Still no paths found, even after successfully adding Retropath edges
        else: # Paths found, return the results
            return paths
    except: # RetroPath failed
        print('retropath failed')
        return False
