import os

from .classes.classes_functions import *
import pandas as pd

def main(start: str, end: str, inchi_start=None, inchi_end=None):
    PATHWAY_FOLDER = os.path.join(os.path.dirname(__file__), 'Pathways')
    os.makedirs(PATHWAY_FOLDER, exist_ok=True)


    # Parameters
    degree_alpha = 2.5  # penalty strength for high-degree nodes
    freq_alpha = 2.5    # penalty strength for frequently-occurring compounds
    MAX_PATHS = 3


    if not os.listdir(PATHWAY_FOLDER):
        raise FileNotFoundError('Warning: the Pathways folder is missing or empty.'
            'Download pathways before attempting to find shortest paths.')

    G_full, G_simple, name_to_ik, ik_to_name, info_dict = build_graphs_with_filters(PATHWAY_FOLDER, degree_alpha=degree_alpha, freq_alpha=freq_alpha)
    #print(info_dict['citrate'])
    #print(info_dict['4.2.1.3'])
    # Normalize keys for prompt (because name_to_ik already uses normalized names in builder)
    start = start.strip().lower()
    end = end.strip().lower()

    start_ik = resolve_name_to_inchi(start, name_to_ik)
    end_ik = resolve_name_to_inchi(end, name_to_ik)

    # 1) Check if start or end is missing in G_simple / G_full
    need_retropath = False
    for g in [G_simple, G_full]:
        if start_ik not in g.nodes or end_ik not in g.nodes:
            need_retropath = True
            break

    if not need_retropath:
        # Try searching G_simple first
        paths = k_shortest_paths(G_simple, start_ik, end_ik, ik_to_name, name_to_ik, K=MAX_PATHS)

        # Fallback to G_full if nothing found
        if not paths:
            # No paths found in filtered network. Trying full network...
            paths = k_shortest_paths(G_full, start_ik, end_ik, ik_to_name, name_to_ik, K=MAX_PATHS)
            if not paths:
                need_retropath = True
            else:
                # Yay paths found, return the results
                return paths,info_dict
        else:
            # Yay paths found, return the results
            return paths,info_dict

    if need_retropath :
        # Fallback to RetroPath
        outdir = "retropath_output/"
        rules_file = "rules.csv"
        returncode, stdout, stderr = rp_single_source(sink=end_ik, rules=rules_file, outdir=outdir, name=start, inchi=start_ik, steps=3, structures=100, dmin=6, dmax=16, timeout=60)

        if returncode == 0:
            # Load Retropath edges
            rp_edges_file = outdir + "rp_edges.tsv"
            rp_df = pd.read_csv(rp_edges_file, sep="\t")
            
            # Add edges to G_simple (or to a new graph copy)
            for _, row in rp_df.iterrows():
                G_simple.add_edge(row['std_inchikey_substrate'], row['std_inchikey_product'], enzyme=row['EC number'])
            
            # Re-run path search
            paths = k_shortest_paths(G_simple, start_ik, end_ik, ik_to_name, name_to_ik, K=MAX_PATHS)

            if not paths:
                # Still no paths found, even after successfully adding Retropath edges
                return False
            else:
                # Yay paths found, return the results
                return paths,info_dict
        else:
            # RetroPath failed or returned no solution
            return False

#print(main('citrate','isocitrate'))
