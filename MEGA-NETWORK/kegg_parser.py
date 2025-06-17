import xml.etree.ElementTree as ET
import csv
import re
import requests
import urllib.parse
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.warning')

enzymes_cache: dict[str, tuple[str, str]] = {}
compounds_cache: dict[str, str] = {}
inchi_cache: dict[str, str] = {}

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
    if not inchi.startswith("InChI="): return ""
    mol = Chem.MolFromInchi(inchi)
    return Chem.MolToInchi(mol, options='-SUCHI') if mol else ""


def get_compound_inchi(name: str) -> str:
    """Resolve a compound name to its standard InChI via PubChem and OPSIN.
    Logs any names for which an InChI could not be found."""
    enc = urllib.parse.quote(name)
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    # 1. Exact PubChem (StandardInChI TXT, then JSON fallback)
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
            std = standardize(inchi)
            if std:
                return std
        except Exception:
            continue

    # 2. OPSIN systematic name parser
    try:
        r = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{enc}.inchi", timeout=8)
        if r.ok:
            std = standardize(r.text.strip())
            if std:
                return std
    except Exception:
        pass

    # 3. PubChem fuzzy lookup (CID list then synonyms)
    try:
        txt = requests.get(f"{base}/name/{enc}/cids/TXT", timeout=8).text
        cids = [cid for cid in txt.split() if cid.isdigit()]
        if cids:
            cid = cids[0]
        else:
            info = requests.get(
                f"{base}/name/{enc}/synonyms/JSON", timeout=8).json()
            tgt = ''.join(ch for ch in name.lower() if ch.isalnum())
            cid = None
            for entry in info.get('InformationList', {}).get('Information', []):
                if any(''.join(ch for ch in s.lower() if ch.isalnum()) == tgt
                       for s in entry.get('Synonym', [])):
                    cid = entry['CID']
                    break
            if cid is None:
                raise ValueError("No CID match via synonyms")

        r2 = requests.get(
            f"{base}/cid/{cid}/property/StandardInChI/TXT", timeout=8)
        r2.raise_for_status()
        std = standardize(r2.text.strip())
        if std:
            return std
    except Exception:
        pass
    return ""


def get_ec_from_gene(gene_ref: str) -> tuple[str, str]:
    """Retrieve enzyme name and EC from an enzyme coding gene's ID"""
    if gene_ref in enzymes_cache:
        return enzymes_cache[gene_ref]
    # Gene → KO
    try:
        r1 = requests.get(f'http://rest.kegg.jp/link/ko/{gene_ref}', timeout=10)
        r1.raise_for_status()
        ko_line = r1.text.strip().split('\n')[0]
        ko_ref = ko_line.split('\t')[1]  # e.g. 'ko:K01080'
    except Exception:
        enzymes_cache[gene_ref] = (gene_ref, '')
        return gene_ref, ''
    # KO → enzyme name and EC
    name = ko_ref
    ec = ''
    try:
        r2 = requests.get(f'http://rest.kegg.jp/get/{ko_ref}', timeout=10)
        r2.raise_for_status()
        for line in r2.text.splitlines():
            if line.startswith('NAME'):
                parts = line.split('NAME', 1)[1].strip()
                name = parts.split('[EC:', 1)[0].strip()
                ec = parts.split('[EC:', 1)[1].split(']', 1)[0].strip()
    except Exception:
        pass
    enzymes_cache[gene_ref] = (name, ec)
    return name, ec


def parse_kgml_to_tsv(kgml_file: str) -> None:
    tree = ET.parse(kgml_file)
    root = tree.getroot()

    # Prepare output filenames
    path_code = root.get('name', '').replace('path:', '')
    nodes_file = f"Pathways/{path_code}_nodes.tsv"
    edges_file = f"Pathways/{path_code}_edges.tsv"

    # Extract <entry> elements of type "compound" or "gene"
    entries: dict[str, dict] = {}
    for e in root.findall('.//{*}entry'):
        etype = e.get('type')
        if etype not in ('compound', 'gene'):
            continue
        eid = e.get('id')
        id2 = e.get('name', '')
        reactions = e.get('reaction', '').split()
        entries[eid] = {'type': etype, 'id2': id2, 'reactions': reactions}

    # Map reactions to enzymes via KO
    reaction2enz: dict[str, list[str]] = {}
    enzyme_info: dict[str, tuple[str, str]] = {}
    missing_ec_file = "missing_ecs.txt"
    with open(missing_ec_file, 'a') as missing_log:
        for eid, info in entries.items():
            if info['type'] == 'gene':
                gene_ref = info['id2'].split()[0]
                enz_name, ec = get_ec_from_gene(gene_ref)

                if ec:  # EC number is present
                    enzyme_info[gene_ref] = (enz_name, ec)
                    for rid in info['reactions']:
                        reaction2enz.setdefault(rid, []).append(gene_ref)
                else:
                    print(f"Skipping gene {gene_ref}: no EC number found")
                    missing_log.write(f"{gene_ref}\tfrom pathway {path_code}\n")


    # Store compounds' names and InChI codes
    compound_entries: dict[str, tuple[str, str]] = {}
    for eid, info in entries.items():
        if info['type'] == 'compound' and info['id2'].startswith('cpd:'):
            c_id = info['id2'].split(':', 1)[1]
            pre_name = get_compound_name(c_id)
            name = re.sub(r'^\[.*?\]\s*', '', pre_name).strip()
            inchi = get_compound_inchi(name)
            clean = inchi.strip().strip('"').splitlines()
            primary_inchi = clean[0] if clean else ""
            compound_entries[eid] = (name, primary_inchi)
            if not primary_inchi:
                with open("missing_inchis.txt", 'a') as f:
                    f.write(f"{name}\tfrom pathway {path_code}\n")

    # Write nodes.tsv
    with open(nodes_file, 'w+', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['name', 'inchi'])
        for name, inchi in compound_entries.values():
            w.writerow([name, inchi])

    # Build edges with reversible handling
    triples: set[tuple[str, str, str, str]] = set()
    for rxn in root.findall('.//{*}reaction'):
        rid = rxn.get('name', '')
        reversible = (rxn.get('type', '').lower() == 'reversible')

        subs = [
            s.get('id') for s in rxn.findall('{*}substrate')
            if s.get('id') in compound_entries
        ]
        prods = [
            p.get('id') for p in rxn.findall('{*}product')
            if p.get('id') in compound_entries
        ]

        for gene_ref in reaction2enz.get(rid, []):
            enz_name, ec = enzyme_info.get(gene_ref, (gene_ref, ''))
            for s in subs:
                s_name, _ = compound_entries[s]
                for p in prods:
                    p_name, _ = compound_entries[p]
                    # forward edge
                    triples.add((enz_name, ec, s_name, p_name))
                    # add reverse edge if reaction is reversible
                    if reversible:
                        triples.add((enz_name, ec, p_name, s_name))

    # Write edges.tsv
    with open(edges_file, 'w+', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['enzyme_name', 'EC', 'source_node_name', 'target_node_name'])
        for enz, ec, s, p in sorted(triples):
            w.writerow([enz, ec, s, p])
