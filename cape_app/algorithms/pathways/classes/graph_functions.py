import os,csv
from itertools import islice
import networkx as nx
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from .chem_utils import *
from .rp_classes import *


def _normalize_name(s: str) -> str:
    """
    Return a stable lowercase normalized string for name matching/lookup.
        Args:
            s: The input string (e.g., a compound name) to normalize.
        Returns:
            A stripped, lowercase version of the input string, or an empty string 
            if the input is not a string or is None.
    """
    return s.strip().lower() if s and isinstance(s, str) else ""


def try_float(x):
    """
    Attempts to convert an input value to a float. If conversion fails, 
    returns the original value.
        Args:
            x: The input value to attempt conversion on.
        Returns:
            The value as a float if successful, otherwise the original value.
    """
    try:
        return float(x)
    except Exception:
        return x


def build_graph_from_tsv(G:nx.MultiDiGraph=nx.MultiDiGraph(), 
                         method:str="KEGG", 
                         rp_outdir:str="./cape_app/pipeline_results/rp_results"):
    """
    Takes a NetworkX MultiDiGraph G (or creates a new empty one) and 
    populates it with edges read from tsv files. 
    The source of the TSV files depends on the 'method' parameter, 
    reading either from KEGG pathway data or RetroPath results.
        Args:
            G: An existing NetworkX MultiDiGraph to add edges to. A new empty 
                MultiDiGraph is created if none is provided.
            method: A string indicating the source of the data. 
                Options are "KEGG" or "RetroPath".
            rp_outdir: The path to the output directory containing RetroPath 
                results, used when method is "RetroPath".
        Returns:
            The populated NetworkX MultiDiGraph G.
    """
    if method=="KEGG":
        pathway_folder = "./cape_app/algorithms/pathways/data/Pathways/"
        suffix = "_edges.tsv"
    elif method=="RetroPath":
        pathway_folder = rp_outdir
        suffix = "weighted_edges.tsv"
    for filename in (f for f in os.listdir(pathway_folder) if f.endswith(suffix)):
        with open(pathway_folder + filename, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                source = row["source_inchi_key"]
                target = row["target_inchi_key"]
                # Define attributes (all except inchikeys) (they might be floats or strings)
                attrs = {k: try_float(v) for k, v in row.items() if k not in ("source_inchi_key", "target_inchi_key")}
                G.add_edge(source, target, **attrs) # add edge, nodes will be added automatically
    return G
    
    
def multigraph_to_digraph(Multi_G: nx.MultiDiGraph) -> nx.DiGraph:
    """
    Convert a NetworkX MultiDiGraph to a DiGraph. If multiple edges exist 
    between the same two nodes, only one edge is kept, and its 'weight' 
    attribute is set to the minimum 'weight' found among the parallel edges 
    in the MultiDiGraph. Node attributes are preserved.
        Args:
            Multi_G: The input NetworkX MultiDiGraph.
        Returns:
            A NetworkX DiGraph with a single edge between any two nodes, 
            representing the minimum weight connection.
    """
    Di_G = nx.DiGraph()
    for u, v, data in Multi_G.edges(data=True):
        weight = float(data.get('weight', 1.0))
        if Di_G.has_edge(u, v):
            Di_G[u][v]['weight'] = min(Di_G[u][v]['weight'], weight)
        else:
            Di_G.add_edge(u, v, weight=weight)
    # Copy node attributes
    for n, attrs in Multi_G.nodes(data=True):
        Di_G.add_node(n, **attrs)
    return Di_G


def k_shortest_paths(G: nx.MultiDiGraph, start_inchikey: str, goal_inchikey: str, K: int = 3):
    """
    Finds the K shortest simple paths between a start and a goal compound 
    (identified by InChIKeys) in a metabolic network represented by a 
    MultiDiGraph. Paths are determined by the minimum edge weight 
    in an intermediate DiGraph representation.
        Args:
            G: The input NetworkX MultiDiGraph representing the metabolic network.
            start_inchikey: The InChIKey of the starting compound.
            goal_inchikey: The InChIKey of the target compound.
            K: The maximum number of shortest simple paths to return. Defaults to 3.
        Returns:
            A list of tuples, where each tuple contains:
            1. The total path weight (float).
            2. A list of path steps. Each step is a list formatted as 
            `[source_name, ec_dicts, target_name]`, where `ec_dicts` is a 
            list of dictionaries containing information about the EC numbers 
            (e.g., {"EC": "1.1.1.1", "enzyme": "...", "color": "...", "reactions": [...]}).
            Paths are sorted by total weight. Returns an empty list if no path is found 
            or if the compounds are not in the graph.
    """
    Di_G = multigraph_to_digraph(G)

    if start_inchikey not in Di_G or goal_inchikey not in Di_G:
        return []

    try:
        path_generator = nx.shortest_simple_paths(Di_G, start_inchikey, goal_inchikey, weight="weight")
    except (nx.NetworkXNoPath, Exception):
        return []

    results = []
    for path in islice(path_generator, K):
        total_w = 0.0
        path_steps = []

        for u, v in zip(path, path[1:]):
            # Weight is calculated from DiGraph (lightest edge)
            w = Di_G[u][v]["weight"]
            total_w += w

            # Collect all parallel edges in the original MultiDiGraph
            edges = G.get_edge_data(u, v, default={})
            ec_dicts_map = {}
            for _, edge_data in edges.items():
                # Get EC info for EC dictionaries
                ec = edge_data.get("EC", "noec")
                enzyme = edge_data.get("enzyme_name", "")
                color = edge_data.get("color", "")
                reaction = edge_data.get("reaction", "")
                reaction_smiles = edge_data.get("Reaction SMILES", "")

                if ec not in ec_dicts_map:
                    ec_dicts_map[ec] = {
                        "EC": ec,
                        "enzyme": enzyme,
                        "color": color,
                        "reactions": [],
                        "Reaction SMILES" : reaction_smiles
                    }
                if reaction:
                    ec_dicts_map[ec]["reactions"].append(reaction)

            # Keep only green if any green exists
            ec_dicts = list(ec_dicts_map.values())
            greens = [d for d in ec_dicts if d["color"] == "green"]
            if greens:
                ec_dicts = greens

            # Use source_name and target_name from any edge because they should be the same
            any_edge = next(iter(edges.values()))
            path_steps.append([any_edge.get("source_name"), ec_dicts, any_edge.get("target_name")])

        results.append((total_w, path_steps))

    results.sort(key=lambda x: x[0])
    return results