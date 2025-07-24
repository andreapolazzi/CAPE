import os
from classes_functions import MetabolicGraph

def get_3_shortest_paths(start:str,end:str):
    # Check that Pathways folder exists or make one
    PATHWAY_FOLDER = 'Pathways'
    os.makedirs(PATHWAY_FOLDER, exist_ok=True)
    # Load all tsv files to networks
    for fname in os.listdir(PATHWAY_FOLDER):
        if not fname.endswith('_edges.tsv'):
            continue
        MetabolicGraph.from_tsv(fname)
    # Create mega network
    MetabolicGraph.mega = MetabolicGraph.megaGraph()

    # Connect start and/or end to central metabolites if they are not in the nodes
    node_list = list(n.lower() for n in MetabolicGraph.mega.nodes)
    if start not in node_list:
        # Connenct start to central metabolites through RetroPath
        pass
    if end not in node_list:
        # Connenct end to central metabolites through RetroPath
        pass

    # Find up to 3 shortest paths
    all_paths = MetabolicGraph.mega.shortest_paths(start, end)

    # Fallback to RetroPath if path can't be found (not yet fully implemented)
    if not all_paths:

        # Connect start to metaboliti principali with RetroPath
        
        # Connect end to metaboliti principali with RetroPath

        # Build graphs from RP results and join them to mega graph
        start_pway_graph = MetabolicGraph.from_tsv(f'.tsv') # da compilare
        end_pway_graph = MetabolicGraph.from_tsv(f'.tsv') # da compilare
        start_pway_graph.join2megaGraph()
        end_pway_graph.join2megaGraph()

        # Retry finding the shortest path
        all_paths = MetabolicGraph.mega.shortest_paths(start, end)
        if not all_paths:
            return "RetroPath 2.0 couldn't find a path."

    # Return text for each candidate path
    WEIGHT_LABEL = {1: "(KEGG green edge)", 10: "(KEGG white edge)", 100: "(RetroPath)"}
    txt = []

    # Sort paths by total distance
    for idx, (dist, path_nodes, path_ecs, path_rxns) in enumerate(sorted(all_paths, key=lambda x: x[0]), start=1):
        header = f"Path #{idx} (total length: {round(dist)})"
        txt.append(header)

        for i in range(len(path_nodes) - 1):
            u, v = path_nodes[i], path_nodes[i+1]
            ec = path_ecs[i]

            # Find the actual edge object
            chosen = next((e for e in u.outgoing if v in e.products and (e.ec == ec or e.name == ec)), None)
            if not chosen:
                chosen = next((e for e in u.outgoing if v in e.products), None)

            w = chosen.weight if chosen else 1
            label = WEIGHT_LABEL.get(w, "")

            step = f"  {i+1}. {u.name} ---[{ec.replace(' ', ' or ')}]---> {v.name} {label}"
            txt[idx-1] += f"\n{step}"

    return txt