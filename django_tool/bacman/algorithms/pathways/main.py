import os
import re
from .classes.classes_functions import MetabolicGraph
from .classes.rp_classes import RP, RpGraph

def main(start: str, end: str, inchi_start=None, inchi_end=None):
    """ Get 3 shortest paths and return structured list output.
        Each path is a list of steps. Each step is:
        [reactant_name (str), [(ec_str, weight_str), ...], product_name (str)]
        non_green_ecs is a list (one element per path) containing EC strings for non-green edges.
    """
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

    # Fallback to RetroPath if path can't be found (not yet fully implemented)
    if not all_paths:
        # Placeholder: build graphs from RP results and join them to mega graph (implement as needed)
        start_pway_graph = MetabolicGraph.from_tsv(f'.tsv')  # da compilare
        start_pway_graph.join2megaGraph()
        end_pway_graph = MetabolicGraph.from_tsv(f'.tsv')  # da compilare
        end_pway_graph.join2megaGraph()

        # Retry finding the shortest path
        all_paths = MetabolicGraph.mega.shortest_paths(start, end)
        if not all_paths:
            return "RetroPath 2.0 couldn't find a path."

    # Mappings for human-friendly weight names
    WEIGHT_COLOR = {1: "green", 10: "white", 100: "RetroPath"}

    paths_list = [] # final data structure for returning paths
    non_green_ecs = []  # list of lists: one per path of the EC strings for non-green edges

    # Sort paths by total distance and build the nested lists
    for i, (total_w, path_nodes, ec_tuples, reaction_ids) in enumerate(sorted(all_paths, key=lambda x: x[0]), start=1): # for each path
        # path_ecs contains lists that can be of one or more EC tuples, because there can be more than one enzyme that can take you from one compound to another
        # the tuples are (EC string, weight of edge)
        path_steps = []
        non_green_ecs.append([])

        for i in range(len(path_nodes) - 1): # for each step in a path
            u,v = path_nodes[i], path_nodes[i+1]
            ecs = ec_tuples[i]

            edge_data = MetabolicGraph.mega.get_edge_data(u, v, default={})
            reactions = edge_data.get("reactions") or []
            if reactions and isinstance(reactions, list) and reactions:
                meta = reactions[0]
                ec = meta.get("ec") or (ec_tuples[i][0] if i < len(ec_tuples) else "-")
                col = meta.get("color")
            else:
                ec = reaction_ecs[i] if i < len(reaction_ecs) else "-"
                col = Fore.RED
            u_name = (ik_to_name.get(u) or used_graph.nodes[u].get("name") or u)
            v_name = (ik_to_name.get(v) or used_graph.nodes[v].get("name") or v)

            # Find the actual edge object ? there's more than one edge, there's one for each EC
            chosen = next((e for e in u.outgoing if v in e.products and (e.ec == ec or e.name == ec)), None)
            if not chosen:
                chosen = next((e for e in u.outgoing if v in e.products), None)

            w = chosen.weight if chosen else 1
            weight_str = WEIGHT_NAME.get(w, "")
            if weight_str != "green":
                non_green_ecs[i-1].append(ec)

            # If ec_item is missing or '', fall back to the chosen edge's ec name
            if ec is None and chosen is not None:
                ec = getattr(chosen, "ec", None) or getattr(chosen, "name", None)

            # Build list of tuples (ec_str, weight_str)
            ecs_with_weights = [(ec_str, weight_str) for ec_str in ecs]

            # If non-green edge, collect its EC(s) into non_green_ecs for this path
            if w != 1:
                non_green_ecs[i-1].extend(ecs)

            # Each step is [reactant_name, [(ec, weight_str), ...], product_name]
            step_item = [str(u.name), ecs_with_weights, str(v.name)]
            path_steps.append(step_item)

        paths_list.append(path_steps)

    return paths_list, non_green_ecs


if __name__ == '__main__':
    #paths_list, non_green_ecs = main(start="Methyl tert-butyl ether", end="acetyl-coa", inchi="InChI=1S/C5H12O/c1-5(2,3)6-4/h1-4H3")
    paths_list, non_green_ecs = main(start="terephthalate", end="succinyl-coa", inchi="InChI=1S/C5H12O/c1-5(2,3)6-4/h1-4H3")

    print(paths_list)
