import os,csv,subprocess
from itertools import islice
import networkx as nx
import pubchempy as pcp
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def chem_distance(substrate_id, product_id, format="inchi"):
    if format == "smiles":
        mol1 = Chem.MolFromSmiles(substrate_id)
        mol2 = Chem.MolFromSmiles(product_id)
    else:  # default InChI
        mol1 = Chem.MolFromInchi(substrate_id)
        mol2 = Chem.MolFromInchi(product_id)
    
    if mol1 is None or mol2 is None:
        return None  # Invalid SMILES/InChI

    # Generate Morgan fingerprints
    fp1 = morgan_gen.GetFingerprint(mol1)
    fp2 = morgan_gen.GetFingerprint(mol2)
    
    # Calculate similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    # Return transformation significance (1 - similarity)
    return 1 - similarity  # Higher = more chemical change


# Map color labels to base numeric weights and viceversa
COLOR_BASE_WEIGHT = {"green": 1.0, "white": 10.0, "RetroPath": 100.0}
WEIGHT_COLOR_BASE = {1.0: "green", 10.0: "white", 100.0: "RetroPath"}


def _normalize_name(s: str) -> str:
    """Return a stable lowercase normalized string for name matching/lookup."""
    return s.strip().lower() if s and isinstance(s, str) else ""


def build_graphs_with_filters(pathway_folder: str):
    """
    Parse *_edges.tsv files and build a graph G with networkx
    Returns (G, name_to_ik, ik_to_name, info_dict).
    """

    # Graph
    G = nx.MultiDiGraph()

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
    for fname in os.listdir(pathway_folder):
        if not fname.endswith("_edges.tsv"):
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
        base_w = float(COLOR_BASE_WEIGHT.get(color, COLOR_BASE_WEIGHT.get("", 100.0)))

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
                s_mass = float(r.get("source_mol_weight") or 0.0)
            except Exception:
                s_mass = None
            try:
                t_mass = float(r.get("target_mol_weight") or 0.0)
            except Exception:
                t_mass = None

            if s_ik:
                reagents.setdefault(s_ik, {"name": s_name or "", "mass": s_mass})
            if t_ik:
                products.setdefault(t_ik, {"name": t_name or "", "mass": t_mass})

            # Store compound info for info_dict
            if s_name and s_ik:
                s_norm = _normalize_name(s_name)
                info_dict.setdefault(_normalize_name(s_norm), {}).update({
                    "inchikey": info_dict.get(s_norm, {}).get("inchikeys", s_ik),
                    "inchi": info_dict.get(s_norm, {}).get("inchis", s_i),
                    "mass": s_mass,
                    "reactions": list(set(info_dict.get(s_norm, {}).get("reactions", []) + [rid]))
                })
            if t_name and t_ik:
                t_norm = _normalize_name(t_name)
                info_dict.setdefault(_normalize_name(t_norm), {}).update({
                    "inchikey": info_dict.get(t_norm, {}).get("inchikeys", t_ik),
                    "inchi": info_dict.get(t_norm, {}).get("inchis", t_i),
                    "mass": t_mass,
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
            mass = v.get("mass")
            if not G.has_node(ik):
                G.add_node(ik, name=name, mass=mass)
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
                G.add_edge(s_ik, t_ik, reaction=rid, ec=ec, color=color, base_weight=base_w, weight=base_w)

    ### Weight Penalties ###

    # Compute chemical distance penalties

    blacklist = [
    # Energy carriers
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump",
    "ctp", "cdp", "cmp", "ppi", "pi", "inorganic phosphate", "orthophosphate", "diphosphate", "reduced fmn",
    # Electron carriers
    "nad+", "nadh", "nadp+", "nadph", "fad", "fadh2", "fmn", "fmnh2","electron-transferring flavoprotein",
    "reduced acceptor", "acceptor", "[reduced nadph---hemoprotein reductase]", "[oxidized nadph---hemoprotein reductase]",
    "reduced ferredoxin", "oxidized ferredoxin", "ferredoxin",
    # Coenzymes and group carriers
    "coenzyme a", "coa", "acetyl-coa", "s-adenosylmethionine", "sam", "thf", "tetrahydrofolate",
    "biotin", "lipoic acid", "pyridoxal phosphate", "plp", "tpp", "thiamine pyrophosphate", "acyl-carrier protein",
    # Gases and small byproducts
    "co2", "carbon dioxide", "h2o", "water", "nh3", "h+", "o2", "oxygen", "h2o2", "ch4", "formate", "formaldehyde",
    # Sugars and nucleotides in transfer reactions
    "udp-glucose", "cdp-choline", "cmp-neu5ac", "prpp", "ribose-5-phosphate", "dump", "dtmp",
    # Miscellaneous
    "glutathione", "gsh", "gssg", "carnitine", "creatine", "phosphocreatine",
    "ubiquinone", "ubiquinol", "cytochrome c", "ferricytochrome c", "fe2+", "fe3+",
    ]

    # Lookup table for chem distances
    chem_cache = {}

    for u, v, key, data in G.edges(keys=True, data=True):
        base = float(data.get("base_weight", data.get("weight", 100.0)))
        u_name = ik_to_name[u].lower()
        v_name = ik_to_name[v].lower()
        if u_name in blacklist or v_name in blacklist:
            data['weight'] = 1e6
            continue
        key_uv = tuple(sorted((u_name, v_name)))
        if key_uv in chem_cache:
            c_dist = chem_cache[key_uv]
        else:
            u_inchi = info_dict[u_name]['inchi']
            v_inchi = info_dict[v_name]['inchi']
            c_dist = chem_distance(u_inchi, v_inchi)
            chem_cache[key_uv] = c_dist
        data['weight'] = base + (c_dist ** 4) * 1000

    return G, name_to_ik, ik_to_name, info_dict


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
        weight = float(data.get('weight', data.get('base_weight', 100.0)))
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
            w = float(edge_data.get('weight', edge_data.get('base_weight', 100.0)))
            total_w += w

            # Determine color from base_weight (not the penalized 'weight')
            base = float(edge_data.get('base_weight', edge_data.get('weight', 100.0)))
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


def get_inchi_from_name(name: str):
    try:
        compounds = pcp.get_compounds(name, 'name')
        if compounds:
            return compounds[0].inchi
        else:
            return None
    except:
        return None


def resolve_name_to_inchi(name, name_to_ik):
    """Try network first, then PubChem if not found"""
    ik = name_to_ik.get(name)
    if ik:
        return ik
    # Try PubChem
    ik = get_inchi_from_name(name)
    if ik:
        # Name not in network, resolved to InChI via PubChem
        return ik
    else:
        # Cannot resolve the name to an InChI and thus cannot run RetroPath
        return False


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