import xml.etree.ElementTree as ET
import csv,re,requests,urllib.parse

enzymes_cache: dict[str, tuple[str, str]] = {}
compounds_cache: dict[str, str] = {}
inchi_cache: dict[str, str] = {}


### Functions used inside the parser

def get_compound_name(c_id: str) -> str:
    if c_id in compounds_cache:
        return compounds_cache[c_id]
    url = f'https://rest.kegg.jp/get/cpd:{c_id}'
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


def get_compound_inchi(name: str) -> str:
    """Fetches the InChI code for a given compound name via the PubChem PUG REST API.
    Handles URL encoding and fallback to JSON parsing for more complex names."""
    # URL-encode the name to handle spaces and special characters
    encoded_name = urllib.parse.quote(name)
    # InChI request from PubChem REST API
    txt_url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded_name}/property/InChI/TXT")
    try:
        r = requests.get(txt_url, timeout=10)
        r.raise_for_status()
        inchi_txt = r.text.strip()
        if inchi_txt.startswith("InChI="):
            return inchi_txt
    except requests.RequestException:
        pass
    # Use JSON to extract InChI for trickier names
    json_url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{encoded_name}/property/InChI,InChIKey/JSON")
    try:
        r2 = requests.get(json_url, timeout=10)
        r2.raise_for_status()
        data = r2.json()
        props = data.get('PropertyTable', {}).get('Properties', [])
        if props:
            inchi = props[0].get('InChI')
            if inchi:
                return inchi
    except (requests.RequestException, ValueError):
        return ""
    return ""


def get_ec_from_gene(gene_ref: str) -> tuple[str, str]:
    """Retrieve enzyme name and EC from an enzyme coding gene's ID"""
    if gene_ref in enzymes_cache:
        return enzymes_cache[gene_ref]
    # Gene → KO
    try:
        r1 = requests.get(f'https://rest.kegg.jp/link/ko/{gene_ref}', timeout=10)
        r1.raise_for_status()
        ko_line = r1.text.strip().split('\n')[0]
        ko_ref = ko_line.split('\t')[1]  # e.g. 'ko:K01080'
    except Exception:
        enzymes_cache[gene_ref] = (gene_ref, '')
        return gene_ref, ''
    # KO → enzyme name and EC
    name = ko_ref
    try:
        r2 = requests.get(f'https://rest.kegg.jp/get/{ko_ref}', timeout=10)
        r2.raise_for_status()
        for line in r2.text.splitlines():
            if line.startswith('NAME'):
                parts = line.split('NAME',1)[1].strip().split(';')
                name = ' '.join(parts[-1].split(' ')[:-1])
                ec = parts[-1].split(' ')[-1][4:-1]
    except Exception:
        pass
    enzymes_cache[gene_ref] = (name, ec)
    return name, ec


 ### Parser for KGML files (from KEGG)

def parse_kgml_to_tsv(kgml_file: str) -> None:
    tree = ET.parse(kgml_file)
    root = tree.getroot()
    # Prepare output file names
    path_code = root.get('name', '').replace('path:', '')
    nodes_file = f"Pathways/{path_code}_nodes.tsv"
    edges_file = f"Pathways/{path_code}_edges.tsv"
    # Extract all entries
    entries: dict[str, dict] = {}
    for e in root.findall('.//{*}entry'):
        etype = e.get('type')
        if etype not in ('compound','gene'):
            continue
        eid = e.get('id')
        id2 = e.get('name','')
        reactions = e.get('reaction','').split()
        entries[eid] = {'type':etype,'id2':id2,'reactions':reactions}
    # Map reactions to enzymes via KO
    reaction2enz: dict[str, list[str]] = {}
    enzyme_info: dict[str, tuple[str,str]] = {}
    for eid,info in entries.items():
        if info['type']=='gene':
            gene_ref = info['id2'].split()[0]
            enz_name, ec = get_ec_from_gene(gene_ref)
            enzyme_info[gene_ref] = (enz_name, ec)
            for rid in info['reactions']:
                reaction2enz.setdefault(rid,[]).append(gene_ref)
    # Store compounds' names and inchi codes
    compound_entries: dict[str, tuple[str,str]] = {}
    for eid,info in entries.items():
        if info['type']=='compound' and info['id2'].startswith('cpd:'):
            c_id = info['id2'].split(':',1)[1]
            pre_name = get_compound_name(c_id) # raw name that can contain the enzyme that the substrate is attached to in [], so we only take the substrate name
            name = re.sub(r'^\[.*?\]\s*', '', pre_name).strip()
            inchi = get_compound_inchi(name)
            try:
                inchi = inchi.split('\n')[0].strip().removeprefix('"').removesuffix('"') # only getting the first inchi because sometimes there is a list
            except:
                pass
            compound_entries[eid] = (name, inchi)
    # Write nodes.tsv
    with open(nodes_file,'w+',newline='',encoding='utf-8') as f:
        w = csv.writer(f,delimiter='\t')
        w.writerow(['name','inchi'])
        for name, inchi in compound_entries.values():
            w.writerow([name, inchi])
    # Build edges
    triples: set[tuple[str,str,str,str]] = set()
    for rxn in root.findall('.//{*}reaction'):
        rid = rxn.get('name','')
        subs = [s.get('id') for s in rxn.findall('{*}substrate') if s.get('id') in compound_entries]
        prods= [p.get('id') for p in rxn.findall('{*}product')   if p.get('id') in compound_entries]
        for gene_ref in reaction2enz.get(rid,[]):
            enz_name, ec = enzyme_info.get(gene_ref,(gene_ref,''))
            for s in subs:
                s_name,_ = compound_entries[s]
                for p in prods:
                    p_name,_ = compound_entries[p]
                    triples.add((enz_name, ec, s_name, p_name))
    # Write edges.tsv
    with open(edges_file,'w+',newline='',encoding='utf-8') as f:
        w = csv.writer(f,delimiter='\t')
        w.writerow(['enzyme_name','EC','source_node_name','target_node_name'])
        for enz, ec, s, p in sorted(triples):
            w.writerow([enz, ec, s, p])