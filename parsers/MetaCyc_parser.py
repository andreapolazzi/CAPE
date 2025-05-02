import xml.etree.ElementTree as ET
import csv

### Parse the file
file_path = "PWY-7794.xml"
tree = ET.parse(file_path)
root = tree.getroot()
ns = {"bp": "http://www.biopax.org/release/biopax-level2.owl#",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#"}

### Find all small molecules, proteins, dna and rna
entity_name = {}
for tag in ("smallMolecule", "protein", "dna", "rna"):
    for elem in root.findall(f".//bp:{tag}", ns):
        entity_id = elem.attrib.get(f"{{{ns['rdf']}}}ID")
        name_elem = elem.find("bp:NAME", ns)
        if entity_id and name_elem is not None: # if the name exists it is stored into the dictionary entity_name
            entity_name[entity_id] = name_elem.text.strip()

### Map every physicalEntityParticipant to an entity ID
pep_to_entity = {}
for pep in root.findall(".//bp:physicalEntityParticipant", ns): # loops through all PEPs
    pep_id = pep.attrib.get(f"{{{ns['rdf']}}}ID")
    phys_ent = pep.find("bp:PHYSICAL-ENTITY", ns)
    if not pep_id or phys_ent is None: # skips if the PEP doesn't have an ID or a corresponding entity
        continue
    ref = phys_ent.get(f"{{{ns['rdf']}}}resource")
    if ref: # if the physical entity is a reference, it saves the ref. ID
        pep_to_entity[pep_id] = ref.lstrip("#")
        continue
    for child in phys_ent.findall("*"): # otherwise it checks for references in child elements
        ref = child.get(f"{{{ns['rdf']}}}resource")
        if ref:
            pep_to_entity[pep_id] = ref.lstrip("#")
            break
    for child in phys_ent: # if there's no ref and the entity is defined inline
        entity_id = child.attrib.get(f"{{{ns['rdf']}}}ID")
        local = child.tag.split("}")[-1]
        if entity_id and local in ("smallMolecule", "protein", "dna", "rna"):
            pep_to_entity[pep_id] = entity_id
            break

### Maps every reaction ID to corresponding enzyme names
reaction_to_enzyme = {}
enzymes = {}
for elem in root.findall(".//bp:catalysis", ns): # loops through all catalysis elements
    eid = elem.attrib.get(f"{{{ns['rdf']}}}ID", "Unknown")
    name_elem = elem.find("bp:NAME", ns)
    if name_elem is not None:
        enzymes[eid] = name_elem.text.strip()
    control = elem.find("bp:CONTROLLED", ns)
    if control is not None:
        target = control.get(f"{{{ns['rdf']}}}resource", "")
        if target.startswith("#"):
            rid = target[1:]
            reaction_to_enzyme[rid] = enzymes.get(eid, eid)

### Function to retrieve substrates or products of a given reaction
def get_reaction_side(reaction, tag):
    names = []
    for pep in reaction.findall(tag, ns): # loops through all elements that are on the left or right side of the reaction (substrates or products)
        ref = pep.get(f"{{{ns['rdf']}}}resource")
        pep_ids = [ref.lstrip("#")] if ref else [
            phys_ent.attrib.get(f"{{{ns['rdf']}}}ID")
            for phys_ent in pep.findall("bp:physicalEntityParticipant", ns)
            if phys_ent.attrib.get(f"{{{ns['rdf']}}}ID")
        ]
        for p_id in pep_ids:
            ent_id = pep_to_entity.get(p_id)
            if not ent_id:
                continue
            name = entity_name.get(ent_id)
            if name:
                names.append(name.strip())
    return names

### Graph: nodes and edges
nodes = set()
edges = []

for reaction in root.findall(".//bp:biochemicalReaction", ns):
    r_id = reaction.attrib.get(f"{{{ns['rdf']}}}ID", "Unknown")
    substrates = get_reaction_side(reaction, "bp:LEFT")
    products = get_reaction_side(reaction, "bp:RIGHT")
    enzyme = reaction_to_enzyme.get(r_id)
    for s in substrates:
        nodes.add((s, "compound"))
        if enzyme and enzyme.strip():
            nodes.add((enzyme, "enzyme"))
            edges.append((s, enzyme, "catalyzed_by"))
    for p in products:
        nodes.add((p, "compound"))
        if enzyme and enzyme.strip():
            edges.append((enzyme, p, "produces"))

# Write nodes.tsv
with open("nodes.tsv", "w", newline='', encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["id", "type"])
    for name, typ in sorted(nodes):
        writer.writerow([name, typ])

# Write edges.tsv
with open("edges.tsv", "w", newline='', encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["source", "target", "interaction"])
    for source, target, interaction in edges:
        writer.writerow([source, target, interaction])