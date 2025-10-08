import re,requests,csv
import xml.etree.ElementTree as ET
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')
from .classes.chem_utils import Std

# Caches
compound_info_cache: dict[str, tuple[str,str]] = {}  # c_id → (name, inchi)
enzymes_cache: dict[str, tuple[str, list[str], list[str]]] = {} # Key: ref_code (KO or gene); Value: (enzyme name, list of ECs, list of reactions)
reaction_dir_cache: dict[str, bool] = {}


def get_compound_info(c_id: str) -> tuple[str, str]:
    """
    Fetches the name and InChI string for a given KEGG compound identifier. 
    Uses a cache to store results. Attempts to retrieve structure information 
    first from KEGG's molfile and falls back to PubChem's SDF via a cross-reference.
        Args:
            c_id: The KEGG compound identifier (e.g., 'C00001').
        Returns:
            A tuple containing the human-readable compound name and the standardized 
            InChI string. If information cannot be found, returns '(c_id, "")'.
    """
    match = re.search(r'C\d{5}', c_id.strip())
    if not match:
        return "Unknown", ""
    c_id = match.group()

    if c_id in compound_info_cache:
        return compound_info_cache[c_id]

    # 1) Fetch KEGG entry to get human-readable name
    try:
        resp = requests.get(f"http://rest.kegg.jp/get/cpd:{c_id}", timeout=10)
        resp.raise_for_status()
    except Exception:
        compound_info_cache[c_id] = (c_id, "")
        return c_id, ""

    name = c_id
    pubchem_cid = None
    for line in resp.text.splitlines():
        if line.startswith("NAME"):
            # Grab first name before any semicolon
            name = line.split("NAME",1)[1].strip().split(";")[0]
        elif line.strip().startswith("PubChem:"):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1].isdigit():
                pubchem_cid = parts[1]

    # 2) Try KEGG molfile → InChI
    inchi = ""
    try:
        mol_resp = requests.get(f"http://rest.kegg.jp/get/cpd:{c_id}/mol", timeout=10)
        mol_resp.raise_for_status()
        mol = Chem.MolFromMolBlock(mol_resp.text, sanitize=True)
        if mol:
            inchi = Chem.MolToInchi(mol, options='-SUCHI')
    except Exception:
        inchi = ""

    # 3) Fallback: PubChem CID → molfile → InChI
    if not inchi and pubchem_cid:
        try:
            sdf_url = (
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
                f"cid/{pubchem_cid}/record/SDF/?record_type=2d"
            )
            sdf_resp = requests.get(sdf_url, timeout=10)
            sdf_resp.raise_for_status()
            mol2 = Chem.MolFromMolBlock(sdf_resp.text, sanitize=True)
            if mol2:
                inchi = Chem.MolToInchi(mol2, options='-SUCHI')
        except Exception:
            inchi = ""

    compound_info_cache[c_id] = (name, inchi)
    return name, inchi


def get_reaction_compounds(rxn_id: str) -> tuple[list[str], list[str]]:
    """
    Retrieves the KEGG compound IDs (C-numbers) for the substrates and 
    products of a given KEGG reaction ID from the KEGG API.
        Args:
            rxn_id: The KEGG reaction identifier (e.g., 'R00001').
        Returns:
            A tuple containing two lists: the list of substrate C-IDs and the 
            list of product C-IDs. Returns '([], [])' on failure.
    """
    url = f"http://rest.kegg.jp/get/{rxn_id}"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
    except Exception as e:
        print(f"Error fetching reaction {rxn_id}: {e}")
        return [], []
    eq = ''
    lines = r.text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('EQUATION'):
            eq = line.split('EQUATION', 1)[1].strip()
            i += 1
            while i < len(lines) and lines[i].startswith(' '):
                eq += ' ' + lines[i].strip()
                i += 1
            break
        i += 1
    if not eq:
        return [], []
    left, right = re.split(r'<=>|=>|<-', eq)
    subs = re.findall(r'C\d{5}', left)
    prods = re.findall(r'C\d{5}', right)
    return subs, prods


# Terms for non-enzymes that don't have an EC

cofactor_terms = ["carrier", "cofactor", "binding protein", "prosthetic group",
    "lipoyl", "biotinyl", "acyl carrier", "thioredoxin", "ferredoxin", "flavodoxin"]

subunit_terms = ["subunit", "beta chain", "alpha chain", "gamma chain", "delta chain",
    "assembly", "scaffold", "component", "domain", "complex"]

regulator_terms = ["regulatory", "activator", "inhibitor", "modulator", "repressor",
    "sensor", "controller", "signal transduction", "stimulator"]

transporter_terms = ["transporter", "channel", "permease", "pump", "port",
    "symporter", "antiporter", "abc transporter", "translocase"]

structural_terms = ["structural", "anchor", "scaffold", "matrix", "cytoskeletal"]


def get_ko_info(ko_ref: str) -> tuple[str, list[str], list[str]]:
    """
    Fetches information (name, EC numbers, and associated KEGG reaction IDs) 
    for a KEGG Orthology (KO) reference. Uses a cache. If no EC number is 
    found, it attempts to infer a non-enzymatic function based on common 
    terms in the name (e.g., 'cofactor', 'transporter').
        Args:
            ko_ref: The KEGG Orthology ID (e.g., 'K00001').
        Returns:
            A tuple containing: 
            1. Enzyme/Protein name (str)
            2. List of EC numbers (list[str]) (or pseudo-ECs if non-enzymatic)
            3. List of associated KEGG reaction IDs (list[str])
    """
    if ko_ref in enzymes_cache:
        return enzymes_cache[ko_ref]
    name, ecs, reactions = ko_ref, [], []
    try:
        r2 = requests.get(f'http://rest.kegg.jp/get/{ko_ref}', timeout=10)
        r2.raise_for_status()
        for line in r2.text.splitlines():
            # EC on the same line as NAME
            if line.startswith('NAME'):
                parts = line.split('NAME', 1)[1]
                name = parts.split('[EC:', 1)[0].strip() if '[EC:' in parts else parts.strip()
                for _, (enz_name, ecs, reactions) in enzymes_cache.items():
                    if enz_name == name:
                        return name, ecs, reactions
                ec_matches = re.findall(r'\[EC:(.*?)\]', line)
                if ec_matches:
                    ecs = ec_matches[0].split()
                    continue
            # EC on the line below of NAME
            elif '[EC:' in line and not ecs:
                ec_matches = re.findall(r'\[EC:(.*?)\]', line)
                if ec_matches:
                    ecs = ec_matches[0].split()
            # get reactions
            elif line.startswith('REACTION'):
                parts = line.split()
                reactions = [p for p in parts[1:] if p.startswith("R") and p[1:].isdigit()]
                break
    except Exception:
        pass
    if not ecs:
        name_lower = name.lower()
        if any(term in name_lower for term in cofactor_terms):
            ecs = [f"{name} (cofactor)"]
        elif any(term in name_lower for term in subunit_terms):
            ecs = [f"{name} (subunit)"]
        elif any(term in name_lower for term in regulator_terms):
            ecs =  [f"{name} (regulator)"]
        elif any(term in name_lower for term in transporter_terms):
            ecs = [f"{name} (transporter)"]
        elif any(term in name_lower for term in structural_terms):
            ecs = [f"{name} (structural)"]
        elif "protein" in name_lower:
            ecs = [f"{name} (accessory protein)"]
        else:
            ecs = [f"{name} (no EC)"]
    enzymes_cache[ko_ref] = (name, ecs, reactions)
    return name, ecs, reactions


def get_gene_info(gene_ref: str) -> tuple[str, list[str], list[str]]:
    """
    Maps a KEGG gene identifier to its corresponding KEGG Orthology (KO) ID, 
    and then retrieves the enzyme information using 'get_ko_info'. Uses a 
    cache to store results.
        Args:
            gene_ref: The KEGG gene identifier (e.g., 'b0001').
        Returns:
            The result of 'get_ko_info(ko_ref)'. Returns '(gene_ref, [], [])' on 
            failure to find the KO link.
    """
    if gene_ref in enzymes_cache:
        return enzymes_cache[gene_ref]
    try:
        r1 = requests.get(f'http://rest.kegg.jp/link/ko/{gene_ref}', timeout=10)
        r1.raise_for_status()
        line = r1.text.strip().split('\n')[0]
        if '\t' not in line or not line:
            enzymes_cache[gene_ref] = (gene_ref, [], [])
            return gene_ref, [], []
        ko_ref = line.split('\t')[1]
    except Exception:
        enzymes_cache[gene_ref] = (gene_ref, [], [])
        return gene_ref, [], []
    return get_ko_info(ko_ref)


def is_reversible_in_kegg(rxn_id: str) -> bool:
    """
    Checks the KEGG reaction definition to determine if it is considered 
    reversible. Uses a cache to store the directionality.
        Args:
            rxn_id: The KEGG reaction identifier (e.g., 'R00001').
        Returns:
            True if the reaction definition contains the reversible symbol ('<=>'), 
            False otherwise.
    """
    if rxn_id in reaction_dir_cache:
        return reaction_dir_cache[rxn_id]
    url = f"http://rest.kegg.jp/get/{rxn_id}"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith("DEFINITION") and "<=>" in line:
                reaction_dir_cache[rxn_id] = True
                return True
    except Exception:
        pass
    reaction_dir_cache[rxn_id] = False
    return False



def parse_kgml_to_tsv(kgml_file: str):
    """
    Parses a KEGG Markup Language (KGML) XML file to extract pathway 
    information and generate two TSV files: one for pathway nodes (compounds) 
    and one for pathway edges (reactions). Edges are categorized as 'green' 
    (direct gene involvement) or 'white' (KO-level association).
        Args:
            kgml_file: The file path to the input KGML XML file.
        Returns:
            None. Writes data to 'Pathways/{path}_nodes.tsv' and 'Pathways/{path}_edges.tsv'.
    """
    tree = ET.parse(kgml_file)
    root = tree.getroot()
    path = root.get('name','').replace('path:','')

    # 1) Reaction directionality map for green edges
    reaction_map: dict[str, str] = {}
    for r in root.findall('.//{*}reaction'):
        rid = r.get('name','').split(':',1)[-1]
        reaction_map[rid] = r.get('type','irreversible')

    nodes_f = f"Pathways/{path}_nodes.tsv"
    edges_f = f"Pathways/{path}_edges.tsv"

    # 2) Collect entries
    entries: dict[str, dict] = {}
    for e in root.findall('.//{*}entry'):
        et, eid = e.get('type'), e.get('id')
        if et in ('compound','gene','ortholog'):
            entries[eid] = {
                'type': et,
                'name': e.get('name',''),
                'reactions': e.get('reaction','').split()
        }

    # 3) Pre‐collect compounds
    compounds: dict[str, tuple[str,str]] = {}
    for info in entries.values():
        if info['type']=='compound' and info['name'].startswith('cpd:'):
            cid = info['name'].split(':',1)[1]
            compounds[cid] = get_compound_info(cid)

    edge_col_values: set[tuple[str,str,str,str,str,str,str,str,str,str]] = set()

    # 4) GREEN edges: genes from KGML reactions
    for eid, info in entries.items():
        if info['type'] != 'gene':
            continue
        pathway_reactions = [rxn.split(':',1)[1] for rxn in info['reactions']]
        for gene_ref in info['name'].split():
            enz_name, ecs, _ = get_gene_info(gene_ref)
            ec = ecs[0] if ecs else ''
            for r_code in pathway_reactions:
                subs, prods = get_reaction_compounds(r_code)
                rev = (reaction_map.get(r_code) == 'reversible')
                for s in subs:
                    s_name, s_inchi = get_compound_info(s)
                    s_inchi_key = Std.standardize(s_inchi)
                    compounds[s] = (s_name, s_inchi)
                    for p in prods:
                        p_name, p_inchi = get_compound_info(p)
                        p_inchi_key = Std.standardize(p_inchi)
                        compounds[p] = (p_name, p_inchi)
                        edge_col_values.add((enz_name, ec, s_inchi, s_inchi_key, p_inchi, p_inchi_key, s_name, p_name, r_code, 'green'))
                        if rev:
                            edge_col_values.add((enz_name, ec, p_inchi, p_inchi_key, s_inchi, s_inchi_key, p_name, s_name, r_code, 'green'))

    # 5) WHITE edges: orthologs from KEGG Reaction definition
    for eid, info in entries.items():
        if info['type']=='ortholog':
            for ko in info['name'].split():
                if ko.startswith('ko:'):
                    name, ecs, reactions = get_ko_info(ko)
                    enz, ec = (name, ecs[0]) if name and ecs else (ko, '')
                    for rid_full in reactions:
                        r_code = rid_full.split(':',1)[-1]
                        subs, prods = get_reaction_compounds(r_code)
                        rev = is_reversible_in_kegg(r_code)
                        for s in subs:
                            s_name, s_inchi = get_compound_info(s)
                            s_inchi_key = Std.standardize(s_inchi)
                            compounds[s] = (s_name, s_inchi)
                            for p in prods:
                                p_name, p_inchi = get_compound_info(p)
                                p_inchi_key = Std.standardize(p_inchi)
                                compounds[p] = (p_name, p_inchi)
                                edge_col_values.add((enz, ec, s_inchi, s_inchi_key, p_inchi, p_inchi_key, s_name, p_name, r_code, 'white'))
                                if rev:
                                    edge_col_values.add((enz, ec, p_inchi, p_inchi_key, s_inchi, s_inchi_key, p_name, s_name, r_code, 'white'))

    # 6) Write nodes.tsv
    with open(nodes_f, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['inchi','name'])
        for cid in sorted(compounds):
            name, inchi = compounds[cid]
            inchi_key = Std.standardize(inchi)
            if not inchi:
                continue
            w.writerow([inchi, inchi_key, name])

    # 7) Write edges.tsv – update header and writerow:
    with open(edges_f, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        # add 'reaction' before 'color'
        w.writerow(['enzyme_name', 'EC', 'source_inchi', 'source_inchi_key', 'target_inchi', 'target_inchi_key', 'source_name', 'target_name', 'reaction', 'color'])
        for en, ec, s_i, s_ik, p_i, p_ik, s_n, p_n, rxn, color in sorted(edge_col_values):
            w.writerow([en, ec, s_i, s_ik, p_i, p_ik, s_n, p_n, rxn, color])