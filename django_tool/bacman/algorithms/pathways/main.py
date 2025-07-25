import os
from .classes.classes_functions import MetabolicGraph
from .classes.rp_classes import RP, RpGraph

def main(start:str, end:str, inchi):
    """ Get 3 shortest paths """
    # Check that Pathways folder exists or make one
    #PATHWAY_FOLDER = '/work/igem/iGEM/django_tool/bacman/algorithms/pathways/Pathways'
    PATHWAY_FOLDER = './bacman/algorithms/pathways/Pathways'  # to run in local
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
    #print(node_list)

    if start not in node_list:
        print("start not in node_list")
        """
        RP.rp_single_source(sink='/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_data/kegg_rp.csv',
                            rules='/work/igem/retrorules_rr02_rp2_flat_forward.csv',
                            outdir='/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results',
                            name=start,
                            inchi=inchi,
                            steps=7,)
        print("RetroPath prediction completed")
        RpGraph.tsv_generator(file_path=f'/work/igem/iGEM/django_tool/bacman/algorithms/pathways/rp_results/Methyl_tert-butyl_ether_scope.csv',
                              outdir='/work/igem/iGEM/django_tool/bacman/algorithms/pathways/Pathways/')
        """
        print("Generating a graph from RetroPath output")
        start_pway_graph = MetabolicGraph.from_tsv('rp_edges.tsv')
        start_pway_graph.join2megaGraph()
        print("Merged the two graphs")

    if end not in node_list:
        print("end not in nodes")
        # DEFINE HOW TO CALL RETROPATH IN THIS CASE
        # end_pway_graph = MetabolicGraph.from_tsv(f'.tsv') # da compilare
        # end_pway_graph.join2megaGraph()


    # Find up to 3 shortest paths
    all_paths = MetabolicGraph.mega.shortest_paths(start, end)
    #print("Finding paths")
    #all_paths = MetabolicGraph.mega.shortest_paths(inchi, end)

    # Fallback to RetroPath if path can't be found (not yet fully implemented)
    if not all_paths:

        # Connect start to metaboliti principali with RetroPath
        
        # Connect end to metaboliti principali with RetroPath

        # Build graphs from RP results and join them to mega graph
        # DEFINE HOW TO CALL RP IN THIS CASE
        start_pway_graph = MetabolicGraph.from_tsv(f'.tsv') # da compilare
        start_pway_graph.join2megaGraph()
        # DEFINE HOW TO CALL RP IN THIS CASE
        end_pway_graph = MetabolicGraph.from_tsv(f'.tsv') # da compilare
        end_pway_graph.join2megaGraph()

        # Retry finding the shortest path
        all_paths = MetabolicGraph.mega.shortest_paths(start, end)
        if not all_paths:
            return "RetroPath 2.0 couldn't find a path."

    # Return text for each candidate path
    WEIGHT_LABEL = {1: "(KEGG green edge)", 10: "(KEGG white edge)", 100: "(RetroPath)"}
    txt = []
    non_green_ecs = []

    # Sort paths by total distance
    for idx, (dist, path_nodes, path_ecs, path_rxns) in enumerate(sorted(all_paths, key=lambda x: x[0]), start=1):
        header = f"Path #{idx} (total length: {round(dist)})"
        txt.append(header)
        non_green_ecs.append([])

        for i in range(len(path_nodes) - 1):
            u, v = path_nodes[i], path_nodes[i+1]
            ec = path_ecs[i]

            # Find the actual edge object
            chosen = next((e for e in u.outgoing if v in e.products and (e.ec == ec or e.name == ec)), None)
            if not chosen:
                chosen = next((e for e in u.outgoing if v in e.products), None)

            w = chosen.weight if chosen else 1
            label = WEIGHT_LABEL.get(w, "")
            if w != 1:
                non_green_ecs[idx-1].append(ec)

            step = f"  {i+1}. {u.name} ---[{ec.replace(' ', ' or ')}]---> {v.name} {label}"
            txt[idx-1] += f"\n{step}"

    return txt, non_green_ecs


if __name__ == '__main__':
    #txt, non_green_ecs = main(start="Methyl tert-butyl ether", end="acetyl-coa", inchi="InChI=1S/C5H12O/c1-5(2,3)6-4/h1-4H3")
    txt, non_green_ecs = main(start="terephthalate", end="succinyl-coa", inchi="InChI=1S/C5H12O/c1-5(2,3)6-4/h1-4H3")
    print(txt)