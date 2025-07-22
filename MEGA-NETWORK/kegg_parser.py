import xml.etree.ElementTree as ET
import csv
import re
import requests
import urllib.parse
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')

# Caches
enzymes_cache: dict[str, tuple[str, list[str], list[str]]] = {} # Key: ref_code (KO or gene); Value: (enzyme name, list of ECs, list of reactions)
enzyme_name_to_ref: dict[str, str] = {}
compounds_cache: dict[str, str] = {}
inchi_cache: dict[str, str] = {}
reaction_dir_cache: dict[str, bool] = {}


def get_compound_name(c_id: str) -> str:
    if c_id in compounds_cache:
        return compounds_cache[c_id]
    url = f'http://rest.kegg.jp/get/cpd:{c_id}'
    name = c_id
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith('NAME'):
                name = line.split('NAME', 1)[1].strip().split(';')[0].strip()
                break
    except Exception:
        pass
    compounds_cache[c_id] = name
    return name

def standardize(inchi: str) -> str:
    if not inchi.startswith("InChI="):
        return ""
    mol = Chem.MolFromInchi(inchi)
    return Chem.MolToInchi(mol, options='-SUCHI') if mol else ""


def get_reaction_compounds(rxn_id: str) -> tuple[list[str], list[str]]:
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
    parts = re.split(r'<=>|=>|<-', eq)
    if len(parts) < 2:
        return [], []
    subs = [comp.strip().split(':')[-1] for comp in parts[0].split('+') if comp.strip()]
    prods = [comp.strip().split(':')[-1] for comp in parts[1].split('+') if comp.strip()]
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


def get_compound_inchi(name: str) -> str:
    enc = urllib.parse.quote(name)
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    # 1. Exact PubChem lookup
    for endpoint in (
        f"{base}/name/{enc}/property/StandardInChI/TXT",
        f"{base}/name/{enc}/property/InChI,InChIKey/JSON",
    ):
        try:
            r = requests.get(endpoint, timeout=8)
            r.raise_for_status()
            data = r.text if endpoint.endswith('.TXT') else r.json()
            inchi = (data.strip() if isinstance(data, str)
                     else data['PropertyTable']['Properties'][0]['InChI'])
            inchi = standardize(inchi)
            if inchi:
                inchi_cache[name] = inchi
                return inchi
        except Exception:
            continue
    # 2. OPSIN
    try:
        r = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{enc}.inchi", timeout=8)
        if r.ok:
            inchi = standardize(r.text.strip())
            if inchi:
                inchi_cache[name] = inchi
                return inchi
    except Exception:
        pass
    # 3. PubChem fuzzy lookup
    try:
        txt = requests.get(f"{base}/name/{enc}/cids/TXT", timeout=8).text
        cids = [cid for cid in txt.split() if cid.isdigit()]
        cid = cids[0] if cids else None
        if cid is None:
            info = requests.get(f"{base}/name/{enc}/synonyms/JSON", timeout=8).json()
            tgt = ''.join(ch for ch in name.lower() if ch.isalnum())
            for entry in info.get('InformationList', {}).get('Information', []):
                if any(''.join(ch for ch in s.lower() if ch.isalnum()) == tgt
                       for s in entry.get('Synonym', [])):
                    cid = entry['CID']; break
        if cid:
            r2 = requests.get(f"{base}/cid/{cid}/property/StandardInChI/TXT", timeout=8)
            r2.raise_for_status()
            inchi = standardize(r2.text.strip())
            if inchi:
                inchi_cache[name] = inchi
                return inchi
    except Exception:
        pass
    return ""


def is_reversible_in_kegg(rxn_id: str) -> bool:
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
    compounds: dict[str, str] = {}
    for eid, info in entries.items():
        if info['type']=='compound' and info['name'].startswith('cpd:'):
            cid = info['name'].split(':',1)[1]
            compounds[eid] = get_compound_name(cid)

    triples: set[tuple[str,str,str,str,str]] = set()

    # 4) GREEN edges: genes from KGML reactions
    for eid, info in entries.items():
        if info['type'] != 'gene':
            continue
        pathway_reactions = [rxn.split(':',1)[1] for rxn in info['reactions']]
        for gene_ref in info['name'].split():
            enz_name, ecs, _ = get_gene_info(gene_ref)
            ec = ecs[0] if ecs else ''
            for code in pathway_reactions:
                subs, prods = get_reaction_compounds(code)
                rev = (reaction_map.get(code) == 'reversible')
                for s in subs:
                    for p in prods:
                        triples.add((enz_name, ec, get_compound_name(s), get_compound_name(p), 'green'))
                        if rev:
                            triples.add((enz_name, ec, get_compound_name(p), get_compound_name(s), 'green'))

    # 5) WHITE edges: orthologs from KEGG Reaction definition
    for eid, info in entries.items():
        if info['type']=='ortholog':
            for ko in info['name'].split():
                if ko.startswith('ko:'):
                    name, ecs, reactions = get_ko_info(ko)
                    enz, ec = (name, ecs[0]) if name and ecs else (ko, '')
                    for rid_full in reactions:
                        code = rid_full.split(':',1)[-1]
                        subs, prods = get_reaction_compounds(code)
                        rev = is_reversible_in_kegg(code)
                        for s in subs:
                            for p in prods:
                                triples.add((enz, ec, get_compound_name(s), get_compound_name(p), 'white'))
                                if rev:
                                    triples.add((enz, ec, get_compound_name(p), get_compound_name(s), 'white'))

    # 6) Write nodes.tsv
    with open(nodes_f, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['name','inchi'])
        for nm in sorted(set(compounds.values())):
            w.writerow([nm, get_compound_inchi(nm)])

    # 7) Write edges.tsv
    with open(edges_f, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['enzyme_name','EC','source','target','color'])
        for en, ec, s, p, color in sorted(triples):
            w.writerow([en, ec, s, p, color])