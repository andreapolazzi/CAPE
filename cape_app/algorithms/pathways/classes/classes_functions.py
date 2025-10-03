import os,csv,subprocess
from itertools import islice
import networkx as nx
from rdkit import DataStructs
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from .chem_utils import *
from .rp_classes import *


# Map color labels to base numeric weights and viceversa
COLOR_BASE_WEIGHT = {"green": 2.0, "white": 4.0, "RetroPath": 16.0}
WEIGHT_COLOR_BASE = {2.0: "green", 4.0: "white", 16.0: "RetroPath"}


def _normalize_name(s: str) -> str:
    """Return a stable lowercase normalized string for name matching/lookup."""
    return s.strip().lower() if s and isinstance(s, str) else ""

def rp_single_source(sink: str, rules: str, outdir: str, name: str, inchi: str, steps=3, structures=100, dmin=6, dmax=16, timeout=60):
    """Launches retropath2 using source name and InChI inside retropath_env 
    Return codes:
        0: No error
        1: File is not found
        2: Running the RetroPath2.0 Knime program produced an OSError
        3: InChI is malformated
        10: Source has been found in the sink (warning)
        11: No solution is found (warning)
        99: unexpected error (not defined in retropath module)"""
    command = f"""
    source /work/igem/anaconda3/etc/profile.d/conda.sh && conda activate retropath_env && \
    python -m retropath2_wrapper {sink} {rules} {outdir} \
    --source_name "{name}" \
    --source_inchi "{inchi}" \
    --max_steps {steps} \
    --topx {structures} \
    --dmin {dmin} \
    --dmax {dmax} \
    --msc_timeout {timeout}
    """
    try:
        result = subprocess.run(command, shell=True, executable="/bin/bash",
                                capture_output=True, text=True, check=False)
        return result.returncode, result.stdout, result.stderr
    except Exception as e:
        return 99, "", str(e)


def try_float(x):
        try:
            return float(x)
        except Exception:
            return x


def build_graph_from_tsv(G:nx.MultiDiGraph=nx.MultiDiGraph(), method:str="KEGG", rp_outdir:str="../rp_results"):
    """
    Takes a nx MultiDiGraph G (or creates a new empty one) and adds edges from the KEGG tsv files or form the RP results (depending on the method param).
    Returns G
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pathway_folder = os.path.join(script_dir, "../data/Pathways")
    for filename in (f for f in os.listdir(pathway_folder) if f.endswith("_edges.tsv")):
        with open(os.path.join(script_dir, "../data/Pathways/")+filename, newline="") as f:
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
    Convert MultiDiGraph to DiGraph keeping only the minimum edge weight.
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
    Returns up to K paths sorted by total weight (using the lightest edges in Di_G).
    Each step is [source, ec_dicts, target]
    where ec_dicts is a list of dicts containing info on ECs, structured as:
      { "EC": "1.1.1.1", "enzyme": "Alcohol dehydrogenase", "color": "green", "reactions": ["R0001", "R0002"] }
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
                        "reaction_SMILES" : reaction_smiles
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

