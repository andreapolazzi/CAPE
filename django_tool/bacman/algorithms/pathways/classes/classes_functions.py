import os,csv,math,subprocess,sys
from itertools import islice
import networkx as nx
import pubchempy as pcp
from colorama import init, Fore, Style


# Cofactors / small molecules
cofactors = {
    "h2o", "co2", "o2", "h+", "pi", "ppi",
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump",
    "nad+", "nadh", "nadp+", "nadph", "fad", "fadh2", "fmn", "fmnh2",
    "coa", "coenzyme a", "sam", "s-adenosyl-l-methionine", "thf", "tetrahydrofolate",
    "nad", "nadp", "nadh2", "fe2+", "fe3+", "h2o2", "nh3", "hco3-",
    "orthophosphate", "inorganic phosphate", "pi", "diphosphate", "pyrophosphate"
}

# Central metabolites
central_metabolites = {
    "glucose", "pyruvate", "phosphoenolpyruvate", "acetyl-coa",
    "succinate", "oxaloacetate", "glyceraldehyde 3-phosphate",
    "fructose-6-phosphate", "dihydroxyacetone phosphate", "phosphoenolpyruvate",
    "d-glyceraldehyde 3-phosphate", "d-glyceraldehyde-3-phosphate",
}

# Map color labels to base numeric weights and viceversa
COLOR_BASE_WEIGHT = {"green": 1.0, "white": 10.0, "RetroPath": 100.0}
WEIGHT_COLOR_BASE = {1.0: "green", 10.0: "white", 100.0: "RetroPath"}


def _normalize_name(s: str) -> str:
    """Return a stable lowercase normalized string for name matching/lookup."""
    return s.strip().lower() if s and isinstance(s, str) else ""


def build_graphs_with_filters(pathway_folder: str,
                              degree_alpha: float = 2.5,
                              freq_alpha: float = 2.5):
    """
    Parse *_edges.tsv files and build two graphs:
      - G_full: connect every reagent -> every product (keeps cofactors in edge metadata)
      - G_simple: simplified connectivity ignoring trivial cofactors (reduces shortcuts)

    Returns (G_full, G_simple, name_to_ik, ik_to_name, info_dict).

    Parameters:
      - cofactors: iterable of names (strings) to treat as trivial cofactors (default DEFAULT_COFACTORS)
      - extra_central_metabolites: iterable of names to treat as central metabolites (gets penalized)
      - degree_alpha / freq_alpha: coefficients controlling penalty strength for degree/frequency
    """
    # Setup cofactor / central metabolite sets

    cofactors_set = set(_normalize_name(x) for x in cofactors)
    central_set = set(_normalize_name(x) for x in central_metabolites)

    # substrings to avoid or penalize strongly (common culprits for weird shortcuts)
    avoid_name_substrings = [
        "strictosidine", "isoipecoside", "puromycin", "alkaloid", "glycoside",
        "ribityl", "lumazine", "flavin", "riboflavin", "ribitol",
        "prpp", "phosphoribosyl", "udp-", "utp-", "dtdp-", "tdp-", "nucleotide", "nucleoside",
        "ribose", "ribulose",
        "-diphospho-", "-pyrophosphate", "-5'-phosphate", "-3'-phosphate",
        "isopren", "terpen", "glycosyl",
        "anthranilate", "diphosphate", "indole"
    ]

    # EC prefixes to penalize (not banned, just discouraged)
    ec_penalty_prefixes = {
        "2.1.1.": 50.0,   # O-/N-methyltransferases
        "2.8.": 30.0,     # sulfur transferases
        "3.2.1.": 40.0,   # glycosidases
    }

    # Graphs to return
    G_full = nx.MultiDiGraph()
    G_simple = nx.MultiDiGraph()

    # Dictionaries for conversion between name and inchikey
    name_to_ik = {}
    ik_to_name = {}

    # Dictionary that collects info about compounds and ECs
    # compound_name: [inchikey, other info that i will add later]
    # EC: [reaction_id, other info that i will add later]
    info_dict = {}
    
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
                    "mass": s_mass,
                    "reactions": list(set(info_dict.get(s_norm, {}).get("reactions", []) + [rid]))
                })
            if t_name and t_ik:
                t_norm = _normalize_name(t_name)
                info_dict.setdefault(_normalize_name(t_norm), {}).update({
                    "inchikey": info_dict.get(t_norm, {}).get("inchikeys", t_ik),
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
        for ik, meta in combined.items():
            name = meta.get("name") or ""
            mass = meta.get("mass")
            if not G_full.has_node(ik):
                G_full.add_node(ik, name=name, mass=mass)
            if not G_simple.has_node(ik):
                G_simple.add_node(ik, name=name, mass=mass)
            if name:
                n = _normalize_name(name)
                if n and n not in name_to_ik:
                    name_to_ik[n] = ik
            ik_to_name.setdefault(ik, name or ik)

        # Identify cofactors
        reagent_list = [(ik, reagents[ik]["name"]) for ik in reagents]
        product_list = [(ik, products[ik]["name"]) for ik in products]
        reag_names = {_normalize_name(n or ""): ik for ik, n in reagent_list}
        prod_names = {_normalize_name(n or ""): ik for ik, n in product_list}

        cofactors_in_rxn = set()
        for n, ik in list(reag_names.items()) + list(prod_names.items()):
            if n in cofactors_set:
                cofactors_in_rxn.add(ik)

        # G_full edges
        for s_ik, _ in reagent_list:
            for t_ik, _ in product_list:
                meta = {
                    "reaction_id": rid,
                    "enzyme": enzyme,
                    "ec": ec,
                    "color": color,
                    "base_weight": base_w,
                    "cofactors": sorted(list(cofactors_in_rxn))
                }
                G_full.add_edge(s_ik, t_ik, reactions=[meta], base_weight=base_w, color=color, weight=base_w)

        # G_simple edges (exclude cofactors)
        main_reag = [ik for ik in reag_names.values() if ik not in cofactors_in_rxn]
        main_prod = [ik for ik in prod_names.values() if ik not in cofactors_in_rxn]
        if not main_reag or not main_prod:
            main_reag = [ik for ik in reag_names.values()]
            main_prod = [ik for ik in prod_names.values()]

        for s_ik in main_reag:
            for t_ik in main_prod:
                meta = {
                    "reaction_id": rid,
                    "enzyme": enzyme,
                    "ec": ec,
                    "color": color,
                    "base_weight": base_w,
                    "cofactors": sorted(list(cofactors_in_rxn))
                }
                G_simple.add_edge(s_ik, t_ik, reactions=[meta], base_weight=base_w, color=color, weight=base_w)

    # Compute degree/frequency penalties for G_simple
    freq = {n: G_simple.in_degree(n) + G_simple.out_degree(n) for n in G_simple.nodes()}
    max_freq = max(freq.values()) if freq else 1
    max_deg = max(dict(G_simple.degree()).values()) if G_simple.number_of_nodes() else 1

    for u, v, key, data in G_simple.edges(keys=True, data=True):
        base = float(data.get("base_weight", data.get("weight", 100.0)))
        deg_u, deg_v = G_simple.degree(u), G_simple.degree(v)
        deg_pen = (math.log1p(deg_u) + math.log1p(deg_v)) / (2 * math.log1p(max_deg)) if max_deg > 0 else 0.0
        freq_pen = (freq.get(u,0)/max_freq + freq.get(v,0)/max_freq)/2.0
        central_pen = ((1 if _normalize_name(G_simple.nodes[u].get("name","")) in central_set else 0)
                       + (1 if _normalize_name(G_simple.nodes[v].get("name","")) in central_set else 0))/2.0
        multiplier = 1.0 + degree_alpha*deg_pen + freq_alpha*freq_pen + 3.0*central_pen
        base_weighted = base * multiplier

        extra_penalty = 0.0
        low_u = _normalize_name(G_simple.nodes[u].get("name","") or "")
        low_v = _normalize_name(G_simple.nodes[v].get("name","") or "")
        if any(sub in low_u for sub in avoid_name_substrings): extra_penalty += 150.0
        if any(sub in low_v for sub in avoid_name_substrings): extra_penalty += 150.0

        ec0 = (data.get("reactions")[0].get("ec","").strip() if data.get("reactions") else "")
        for pref, p in ec_penalty_prefixes.items():
            if ec0.startswith(pref):
                extra_penalty += p
                break

        # Mass penalty
        u_mass_val = float(G_simple.nodes[u].get("mass") or 0.0)
        v_mass_val = float(G_simple.nodes[v].get("mass") or 0.0)
        if u_mass_val>1.0 and v_mass_val>0:
            ratio = v_mass_val / u_mass_val
            if ratio < 0.8:
                extra_penalty += (1.0 - ratio) * 120.0
            elif ratio > 1.5:
                extra_penalty += (ratio - 1.0) * 80.0

        # Rarity penalty
        rarity_penalty = 0.0
        if deg_u < 2: rarity_penalty += 20*(2-deg_u)
        if deg_v < 2: rarity_penalty += 20*(2-deg_v)
        if freq.get(u,0)<2: rarity_penalty += 15*(2-freq.get(u,0))
        if freq.get(v,0)<2: rarity_penalty += 15*(2-freq.get(v,0))
        if 0<u_mass_val<50 or u_mass_val>1000: rarity_penalty += 10
        if 0<v_mass_val<50 or v_mass_val>1000: rarity_penalty += 10
        extra_penalty += rarity_penalty

        final_weight = base_weighted + extra_penalty
        data['weight'] = final_weight

    # Ensure G_full edges have weight set
    for u,v,key,data in G_full.edges(keys=True,data=True):
        if "weight" not in data:
            data["weight"] = float(data.get("base_weight",100.0))
        data.setdefault("color", data.get("color",""))
        data.setdefault("reactions", data.get("reactions",[]))

    return G_full, G_simple, name_to_ik, ik_to_name, info_dict


def multigraph_to_digraph(Multi_G: nx.MultiDiGraph) -> nx.DiGraph:
    """
    Convert a MultiDiGraph to a DiGraph by keeping the minimum weight edge between nodes
    and aggregating all reactions into a list.
    """
    Di_G = nx.DiGraph()
    for u, v, data in Multi_G.edges(data=True):
        if Di_G.has_edge(u, v):
            # Keep minimum weight
            Di_G[u][v]['weight'] = min(Di_G[u][v]['weight'], data.get('weight', 100.0))
            # Aggregate reactions
            Di_G[u][v]['reactions'].extend(data.get('reactions', []))
        else:
            Di_G.add_edge(u, v, weight=data.get('weight', 100.0), reactions=list(data.get('reactions', [])))
    # Copy node attributes
    for n, attrs in Multi_G.nodes(data=True):
        Di_G.add_node(n, **attrs)
    return Di_G


def k_shortest_paths(G: nx.DiGraph, start_inchikey: str, goal_inchikey: str, ik_to_name:dict, name_to_ik:dict, K=3, weight="weight"):
    """
    It returns a list of up to K paths, sorted by total weight (smallest first).
    Each element is a tuple: (total_weight:float, path_steps:list).
    Each step in the path is a list: [compound1:str, ec_tuples:list, compound2:str],
    where ec_tuples is a list of tuples (EC:str, color:str, reaction_id:str),
    which represents all enzymatic reactions connecting the two compounds.
    If any green edges exist, only the green ones are included; otherwise, all ECs are included (yellow and red).
    """
    G = multigraph_to_digraph(G)

    if start_inchikey not in G or goal_inchikey not in G:
        return []

    try:
        path_generator = nx.shortest_simple_paths(G, start_inchikey, goal_inchikey, weight=None)
    except nx.NetworkXNoPath:
        # No path exists between nodes"
        return []
    except Exception as e:
        # Path generator failed
        return []

    results = []
    for path in islice(path_generator, K): # For each path
        total_w = 0.0
        path_steps = []
        for u, v in zip(path, path[1:]): # For each step in the current path
            ec_tuples = [] # [(EC, color, reaction_id), ...] # info on all ECs for a given step
            edge_data = G.get_edge_data(u, v, default={})
            edge_weights = []
            if edge_data:  # only if edge exists
                w = edge_data.get("weight", edge_data.get("base_weight", 100.0))
                edge_weights.append(w)
                reactions = edge_data.get("reactions", [])
                for rxn in reactions:
                    ec_tuples.append((rxn.get("ec") or rxn.get("enzyme") or "",
                    WEIGHT_COLOR_BASE.get(rxn.get("base_weight", 100.0), "red"),
                    rxn.get("reaction_id") or ""))

            # For calculating the weight of the step, consider the lightest edge
            if edge_weights:
                total_w += min(edge_weights)

            # For each step, keep only green edges if present
            greens = [t for t in ec_tuples if t[1] == "green"]
            if greens:
                ec_tuples = greens

            # Add this step to the path
            path_steps.append([ik_to_name[u], ec_tuples, ik_to_name[v]])
        results.append((total_w, path_steps))

    # Sort by total weight so that shorter / lighter paths are first
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
        print(f"'{name}' not in network, resolved to InChI via PubChem.")
        return ik
    else:
        print(Fore.RED + f"Cannot resolve '{name}' to an InChI. Cannot run RetroPath." + Style.RESET_ALL)
        sys.exit(1)


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