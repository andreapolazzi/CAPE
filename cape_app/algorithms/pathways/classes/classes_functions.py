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


# questa la ho lasciata come prima perchè non ho capito bene da dove prendere la rp_score
def add_weight_penalties(G:nx.MultiDiGraph, info_dict:dict, ik_to_name:dict):
    blacklist = [
        # Energy carriers / Electron carriers
        "ZKHQWZAMYRWXGA-KQYNXXCUSA-N", # ATP
        "XTWYTFMLZFPYCI-KQYNXXCUSA-N", # ADP
        "UDMBCSSLTHHNCD-KQYNXXCUSA-N", # AMP
        "XKMLYUALXHKNFT-UUOKFMHZSA-N", # GTP
        "QGWNDRXFNXRZMB-UUOKFMHZSA-N", # GDP
        "RQFCJASXJCIDSX-UUOKFMHZSA-N", # GMP
        "PGAVKCOVUIYSFO-XVFCMESISA-N", # UTP
        "XCCTYIAWTASOJW-XVFCMESISA-N", # UDP
        "DJJCXFVJDGTHFX-XVFCMESISA-N", # UMP
        "PCDQPRRSZKQHHS-XVFCMESISA-N", # CTP
        "ZWIADYZPOWUWEW-XVFCMESISA-N", # CDP
        "IERHLVCPSMICTF-XVFCMESISA-N", # CMP
        "NBIIXXVUZAFLBC-UHFFFAOYSA-N", # Orthophosphate
        "XPPKVPWEQAFLFU-UHFFFAOYSA-N", # Diphosphate
        "UNXRWKVEANCORM-UHFFFAOYSA-N", # Polyphosphate
        "YTNIXZGTHTVJBW-SCRDCRAPSA-N", # Reduced FMN
        "BAWFJGJZGIEFAR-NNYOXOHSSA-O", # NAD+
        "BOPGDPNILDQYTO-NNYOXOHSSA-N", # NADH
        "XJLXINKUBYWONI-NNYOXOHSSA-O", # NADP+
        "ACFIXJIJDZMPPO-NNYOXOHSSA-N", # NADPH
        "VWWQXMAJTJZDQX-UYBVJOGSSA-N", # FAD
        "YPZRHBJKEMOYQH-UYBVJOGSSA-N", # FADH2
        "FVTCRASFADXXNN-SCRDCRAPSA-N", # FMN
        # Coenzymes and group carriers
        "RGJOEKWQDUBAIZ-IBOSZNHHSA-N", # CoA, Coenzyme A
        "ZSLZBFCDCINBPY-KMYLAXNMSA-N", # Acetyl-CoA
        "MSTNYGQPCMXVAQ-PXYINDEMSA-N", # Tetrahydrofolate
        "YBJHBAHKTGYVGT-ZKWXMUAHSA-N", # Biotin
        "NGVDGCNFYWLIFO-UHFFFAOYSA-N", # Pyridoxal phosphate
        "AYEKOFBPNLCAJY-UHFFFAOYSA-O", # Thiamine Diphosphate
        # Gases and small byproducts
        "CURLTUGMZLYLDI-UHFFFAOYSA-N", # CO2, Carbon Dioxide
        "XLYOFNOQVPJJNP-UHFFFAOYSA-N", # H2O, Water
        "QGZKDVFQNNGYKY-UHFFFAOYSA-N", # NH3, Ammonia
        "GPRLSGONYQIRFK-UHFFFAOYSA-N", # H+, Hydrogen
        "MYMOFIZGZYHOMD-UHFFFAOYSA-N", # O2, Oxygen
        "XCCTYIAWTASOJW-XVFCMESISA-N", # H2O2, Hydrogen Peroxide
        "VNWKTOKETHGBQD-UHFFFAOYSA-N", # Methane, CH4
        "BDAGIHXWWSANSR-UHFFFAOYSA-N", # Formate
        "WSFSSNUMVMOOMR-UHFFFAOYSA-N", # Formaldehyde
        # Sugars and nucleotides in transfer reactions
        "HSCJRCZFDFQWRP-JZMIEXBBSA-N", # UDP-Glucose
        "RZZPDXZPRHQOCG-OJAKKHQRSA-N", # CDP-Choline
        "JSRLJPSBLDHEIO-SHYZEUOFSA-N", # dUMP
        "GYOZYWVXFNDGLU-XLPZGREQSA-N", # dTMP
        # Miscellaneous
        "RWSXRVCMGQZWBV-UHFFFAOYSA-N", # Glutathione
        "PHIQHXFUZVPYII-UHFFFAOYSA-O", # Carnitine
        "CVSVTCORWBXHQV-UHFFFAOYSA-N", # Creatine
        "DRBBFCLWYRJSJZ-UHFFFAOYSA-N", # Phosphocreatine
        "CRADQHJNYKKORV-UHFFFAOYSA-N", # Ubiquinone
        "RJUKAZDREBHPDA-UHFFFAOYSA-N", # Ubiquinol
        "GAOZTHIDHYLHMS-UHFFFAOYSA-N", # Ferricytochrome b5
        "LJZLHDYZFURWEE-UHFFFAOYSA-N", # Ferricytochrome cL
    ]

    # Lookup table for chemical distances
    dist_cache = {}

    for u, v, key, data in G.edges(keys=True, data=True):
        base = float(data.get("base_weight", data.get("weight", 16.0)))
        u_name = ik_to_name[u].lower()
        v_name = ik_to_name[v].lower()
        if u in blacklist or v in blacklist:
            data['weight'] = 1e6
            continue
        key_uv = tuple(sorted((u_name, v_name)))
        if key_uv in dist_cache:
            dist = dist_cache[key_uv]
        else:
            u_inchi = info_dict[u_name]['inchi']
            v_inchi = info_dict[v_name]['inchi']
            dist = Penalty.transf_distance(u_inchi, v_inchi)
            dist_cache[key_uv] = dist
            d_penalty = Penalty.dist_penalty(dist)
            s_mw,t_mw = data['source_mol_weight'],data['target_mol_weight']
            if s_mw != 0.0 and t_mw != 0.0:
                mw_penalty = Penalty.mw_change(s_mw,t_mw)
            else:
                mw_penalty = 0 # if we don't know the weight idk if we should put a higher difference or 0 ???
        rp_score = data['rp_score']
        final_penalty = Penalty.transf_penalty(rp_score,d_penalty,mw_penalty)
        data['weight'] = final_penalty
    return G


def build_graph_from_tsv(G:nx.MultiDiGraph=nx.MultiDiGraph(), method:str="KEGG"):
    """
    Parse *_edges.tsv files or RetroPath results and build a graph G with networkx or add edges to an existing graph
    Returns (G, name_to_ik, ik_to_name, info_dict).
    """
    
    # Dictionaries for conversion between name and inchikey
    name_to_ik = {}
    ik_to_name = {}

    # Dictionary that collects info about compounds and ECs
    # compound_name: [inchikey, other info that i will add later]
    # EC: [reaction_id, other info that i will add later]
    info_dict = {}

    ### Read tsv files and create networkx graph and fill the dictionaries ###
    
    # Read tsv files and group rows by reaction id
    reactions_store = {}
    if method=="KEGG":
        pathway_folder = "Pathways"
        file_end = "_edges.tsv"
    else: # RetroPath
        pathway_folder = "rp_results"
        file_end = "rp_edges.tsv"
    for fname in os.listdir(pathway_folder):
        if not fname.endswith(file_end):
            continue
        filepath = os.path.join(pathway_folder, fname)
        try:
            with open(filepath, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    rid = (row.get("reaction") or row.get("reaction_id") or "").strip()
                    if not rid:
                        continue
                    reactions_store.setdefault(rid, []).append((filepath, row))
        except Exception:
            continue

    for rid, rows in reactions_store.items():
        _, first = rows[0]
        enzyme = (first.get("enzyme_name") or "").strip()
        ec = (first.get("EC") or "").strip()
        color = (first.get("color") or "").strip().lower()
        base_w = float(COLOR_BASE_WEIGHT.get(color, COLOR_BASE_WEIGHT.get("", 16.0)))

        reagents = {}
        products = {}
        for path, r in rows:
            s_ik = (r.get("source_inchi_key") or "").strip()
            t_ik = (r.get("target_inchi_key") or "").strip()
            s_i = (r.get("source_inchi") or "").strip()
            t_i = (r.get("target_inchi") or "").strip()
            s_name = (r.get("source_name") or "").strip()
            t_name = (r.get("target_name") or "").strip()
            try:
                s_mw = float(r.get("source_mol_weight") or 0.0)
            except Exception:
                s_mw = 0.0
            try:
                t_mw = float(r.get("target_mol_weight") or 0.0)
            except Exception:
                t_mw = 0.0
            try:
                rp_score = float((r.get("RP_score") or 0.0).strip())
            except Exception:
                rp_score = 0.0

            if s_ik:
                reagents.setdefault(s_ik, {"name": s_name or "", "mol_weight": s_mw})
            if t_ik:
                products.setdefault(t_ik, {"name": t_name or "", "mol_weight": t_mw})

            # Store compound info for info_dict
            if s_name and s_ik:
                s_norm = _normalize_name(s_name)
                info_dict.setdefault(_normalize_name(s_norm), {}).update({
                    "inchikey": info_dict.get(s_norm, {}).get("inchikey", s_ik),
                    "inchi": info_dict.get(s_norm, {}).get("inchi", s_i),
                    "mol_weight": s_mw,
                    "reactions": list(set(info_dict.get(s_norm, {}).get("reactions", []) + [rid]))
                })
            if t_name and t_ik:
                t_norm = _normalize_name(t_name)
                info_dict.setdefault(_normalize_name(t_norm), {}).update({
                    "inchikey": info_dict.get(t_norm, {}).get("inchikey", t_ik),
                    "inchi": info_dict.get(t_norm, {}).get("inchi", t_i),
                    "mol_weight": t_mw,
                    "reactions": list(set(info_dict.get(t_norm, {}).get("reactions", []) + [rid]))
                })

        # Also store EC info for info_dict
        if ec:
            info_dict.setdefault(ec, {}).update({
                "enzyme": enzyme,
                "reactions": list(set(info_dict.get(ec, {}).get("reactions", []) + [rid]))
            })

        combined = {**reagents, **products}
        for ik, v in combined.items():
            name = v.get("name") or ""
            mol_weight = v.get("mol_weight")
            if not G.has_node(ik):
                G.add_node(ik, name=name, mol_weight=mol_weight)
            if name:
                n = _normalize_name(name)
                if n and n not in name_to_ik:
                    name_to_ik[n] = ik
            ik_to_name.setdefault(ik, name or ik)

        # Identify cofactors
        reagent_list = [(ik, reagents[ik]["name"]) for ik in reagents]
        product_list = [(ik, products[ik]["name"]) for ik in products]

        # G edges
        for s_ik, _ in reagent_list:
            for t_ik, _ in product_list:
                G.add_edge(s_ik, t_ik, reaction=rid, ec=ec, color=color, source_mol_weight=s_mw, target_mol_weight=t_mw, base_weight=base_w, weight=base_w, rp_score=rp_score)

    return add_weight_penalties(G, info_dict, ik_to_name), name_to_ik, ik_to_name, info_dict


def multigraph_to_digraph(Multi_G: nx.MultiDiGraph) -> nx.DiGraph:
    """
    Convert a MultiDiGraph to a DiGraph while preserving useful metadata:
    - weight (minimum of multi-edges)
    - base_weight (minimum base_weight)
    - ecs (list of ECs / enzyme identifiers)
    - reactions (list of reaction ids)
    - colors (list of color strings from edges)
    """
    Di_G = nx.DiGraph()
    for u, v, data in Multi_G.edges(data=True):
        # Gather fields from multi-edge
        weight = float(data.get('weight', data.get('base_weight', 16.0)))
        base = float(data.get('base_weight', data.get('weight', weight)))
        # Support several possible keys used in different parts of pipeline
        ec = data.get('ec') or data.get('EC') or data.get('enzyme') or ""
        reaction = data.get('reaction') or data.get('reactions') or []
        if isinstance(reaction, str):
            reaction = [reaction] if reaction else []
        color = data.get('color')

        if Di_G.has_edge(u, v):
            # keep minimal numeric weight
            Di_G[u][v]['weight'] = min(Di_G[u][v].get('weight', float('inf')), weight)
            # aggregate reactions
            Di_G[u][v].setdefault('reactions', [])
            for r in reaction:
                if r and r not in Di_G[u][v]['reactions']:
                    Di_G[u][v]['reactions'].append(r)
            # aggregate ecs
            Di_G[u][v].setdefault('ecs', [])
            if ec and ec not in Di_G[u][v]['ecs']:
                Di_G[u][v]['ecs'].append(ec)
            # keep minimal base_weight
            Di_G[u][v]['base_weight'] = min(Di_G[u][v].get('base_weight', base), base)
            # aggregate colors
            Di_G[u][v].setdefault('colors', [])
            if color and color not in Di_G[u][v]['colors']:
                Di_G[u][v]['colors'].append(color)
        else:
            Di_G.add_edge(
                u, v,
                weight=weight,
                base_weight=base,
                ecs=[ec] if ec else [],
                reactions=list(reaction),
                colors=[color] if color else []
            )

    # copy node attributes
    for n, attrs in Multi_G.nodes(data=True):
        Di_G.add_node(n, **attrs)
    return Di_G


def k_shortest_paths(G: nx.DiGraph, start_inchikey: str, goal_inchikey: str, ik_to_name: dict, K: int = 3):
    """
    Returns up to K paths sorted by total weight.
    Each step is [compound_from, ec_tuples, compound_to]
    where ec_tuples is a list of (EC_str, color_str).
    Uses base_weight to decide color (so penalized 'weight' does not change color).
    """
    G = multigraph_to_digraph(G)

    if start_inchikey not in G or goal_inchikey not in G:
        return []

    try:
        # Use weighted shortest_simple_paths by weight
        path_generator = nx.shortest_simple_paths(G, start_inchikey, goal_inchikey, weight='weight')
    except nx.NetworkXNoPath:
        return []
    except Exception:
        return []

    results = []
    for path in islice(path_generator, K):
        total_w = 0.0
        path_steps = []
        for u, v in zip(path, path[1:]):
            edge_data = G.get_edge_data(u, v, default={})
            # Use the actual numeric weight for path cost (already computed in DiGraph)
            w = float(edge_data.get('weight', edge_data.get('base_weight', 16.0)))
            total_w += w

            # Determine color from base_weight (not the penalized 'weight')
            base = float(edge_data.get('base_weight', edge_data.get('weight', 16.0)))
            color_name = WEIGHT_COLOR_BASE.get(base, "red")

            ecs = edge_data.get('ecs', []) or []

            if not ecs:
                ec_tuples = [("", color_name)]
            else:
                ec_tuples = [(ec, color_name) for ec in ecs]

            # If any green entries exist, keep only the green ones
            greens = [t for t in ec_tuples if t[1] == "green"]
            if greens:
                ec_tuples = greens

            path_steps.append([ik_to_name.get(u, u), ec_tuples, ik_to_name.get(v, v)])
        results.append((total_w, path_steps))

    results.sort(key=lambda x: x[0])
    return results
