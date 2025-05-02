import xml.etree.ElementTree as ET
import csv
import requests
import re


### Functions
compounds = {} # to store names of compund IDs that we've already seen
enzymes = {} # to store names of enzyme IDs that we've already seen

def get_compound_name(c_id):
    '''For a given KEGG Compound ID it gets its NAME from KEGG.'''
    if c_id in compounds: # if we've already seen the ID, get the name from the compound dictionary
        return compounds[c_id]
    url = f'https://rest.kegg.jp/get/cpd:{c_id}' # otherwise use the API
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith('NAME'):
                name = line.split('NAME',1)[1].strip().split(';')[0].strip() # gets the first name (sometimes there's more than one)
                break
    except:
        name = c_id # if it hasn't found a name, it uses the compound ID
        pass
    compounds[c_id] = name
    return name

def get_enzyme_name(enz_id):
    '''For a given KEGG Enzyme/Gene ID it gets its NAME from KEGG and it remove the source '(...)'.'''
    if enz_id in enzymes: # if we've already seen the ID, get the name from the enzyme dictionary
        return enzymes[enz_id]
    url = f'https://rest.kegg.jp/get/{enz_id}' # otherwise use the API
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith('NAME'):
                parts = line.split('NAME',1)[1].strip().split(';')
                if len(parts) >= 2: # when there are multiple names, the first is often a short abbreviation
                    name = parts[1].strip() # which is why we pick the second one, which is the 'real' name
                else: # if it's just one name then we just use that one
                    name = parts[0].strip()
                name = re.sub(r'^\([^)]*\)\s*', '', name) # to remove the source, such as '(RefSeq)' or others that can be found before the enzyme name
                break
    except:
        name = enz_id # if it hasn't found a name, it uses the enzyme ID
        pass
    enzymes[enz_id] = name
    return name


### Parse the KGML file to get an ElementTree object
tree = ET.parse('TAG.kgml')
root = tree.getroot()


### Get compounds & genes
entries = {}
for e in root.findall('.//{*}entry'): # loops through every <entry> element
    entry_id   = e.get('id')
    entry_type = e.get('type')
    if entry_type not in ('compound','gene'): # we only want compund or gene (enzyme-coding), because enzyme only has the EC number
        continue
    id2   = e.get('name','') # the compound or gene id
    gfx   = e.find('{*}graphics')
    label = gfx.get('name') if gfx is not None else id2 # looks for the KEGG display name
    reactions = e.get('reaction','').split() # which reactions it's in

    entries[entry_id] = {
        'type':      entry_type,
        'id2':       id2,
        'label':     label,
        'name':      None, # placeholder
        'reactions': reactions
    }


### Get Names
for entry_id,info in entries.items():
    if info['type'] == 'compound' and info['id2'].startswith('cpd:'):
        comp_id = info['id2'].split(':',1)[1]
        info['name'] = get_compound_name(comp_id)
    elif info['type'] == 'gene':
        gene_id = info['id2'].split()[0]
        info['name'] = get_enzyme_name(gene_id)


### Dictionary that maps reaction IDs to gene/enzyme IDs
reaction_to_genes = {}
for eid,info in entries.items():
    if info['type'] == 'gene':
        for reaction in info['reactions']:
            if reaction in reaction_to_genes.keys():
                reaction_to_genes[reaction].append(eid)
            else:
                reaction_to_genes[reaction] = [eid,]


### Edges: substrate to enzyme and enzyme to product, write them in a edges.tsv tab-delim file
edges = []
for reaction in root.findall('.//{*}reaction'): # loops over all <reaction> elements
    r_codes = reaction.get('name','').split()
    genes  = [g for r_code in r_codes for g in reaction_to_genes.get(r_code, [])] # all genes that catalyze the reaction the reaction_to_genes dictionary
    substrates  = [s.get('id') for s in reaction.findall('{*}substrate') if s.get('id') in entries] # get all substrate elements in the reaction, which are in the entries (compound and gene)
    products = [p.get('id') for p in reaction.findall('{*}product')   if p.get('id') in entries] # get all product elements in the reaction, which are in the entries (compound and gene)
    for s in substrates:
        for g in genes:
            edges.append((entries[s]['name'], entries[g]['name'], 'catalyzed_by')) # substrate to enzyme
    for g in genes:
        for p in products:
            edges.append((entries[g]['name'], entries[p]['name'], 'produces')) # enzyme to product

with open('edges.tsv','w',newline='',encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['source','target','interaction'])
    for src, tgt, iact in edges:
        writer.writerow([src, tgt, iact])


### Nodes, write them in a nodes.tsv tab-delim file
nodes = {}
for info in entries.values():
    if info['name'] is None:
        continue
    ntype = 'enzyme' if info['type']=='gene' else 'compound'
    color = 'purple' if ntype=='enzyme' else 'lightblue' # enzymes are colored in purple and compunds are colored in lightblue
    nodes[info['name']] = (ntype, color)

with open('nodes.tsv','w',newline='',encoding='utf-8') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['id','type','color'])
    for name,(ntype,color) in sorted(nodes.items()):
        writer.writerow([name, ntype, color])