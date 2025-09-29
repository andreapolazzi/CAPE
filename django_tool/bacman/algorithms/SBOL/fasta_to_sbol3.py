#!/usr/bin/env python3
"""
Convert one or more FASTA files of parts into a single SBOL3 document.

- Auto-infers Component.type and roles from header keywords and filenames.
- Respects explicit header tags: type=..., role=SO:...
- Creates a Collection to group all components.
- Outputs SBOL3 XML and JSON.

Dependencies: pip install sbol3
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

try:
    import sbol3
except ImportError:
    sys.stderr.write("Missing dependency: pip install sbol3\n")
    sys.exit(1)

# ---------- Heuristics & Ontology helpers ----------

SO = {
    "promoter": "SO:0000167",
    "rbs": "SO:0000139",
    "cds": "SO:0000316",
    "terminator": "SO:0000141",
    "origin_of_replication": "SO:0000296",
    "vector_backbone": "SO:0000755",  # general backbone region
    "operator": "SO:0000057",
    "utr": "SO:0000204",
    "linker": "SO:0001221",
    "scar": "SO:0001953",
}

KEYWORD_TO_ROLE = [
    (r"\b(promoter|pJ231|proD|t7prom)\b", SO["promoter"]),
    (r"\b(rbs|b003[0-9]|shine[- ]?dalgarno)\b", SO["rbs"]),
    (r"\b(cds|gene|orf|coding|cds_|_cds|_gene)\b", SO["cds"]),
    (r"\b(term|terminator|rrnB_T1|T7term)\b", SO["terminator"]),
    (r"\b(backbone|vector|plasmid)\b", SO["vector_backbone"]),
    (r"\b(ori|origin)\b", SO["origin_of_replication"]),
    (r"\b(operator|lacO|tetO)\b", SO["operator"]),
    (r"\b(utr|5utr|3utr)\b", SO["utr"]),
    (r"\b(linker|spacer|gs|gly-ser|gggs)\b", SO["linker"]),
]

def guess_component_type(seq: str, header: str, filename: str) -> str:
    """Return SBOL3 component type URI (DNA/RNA/Protein)."""
    # Allow explicit override in header: type=protein | type=rna | type=dna
    m = re.search(r"\btype=(protein|rna|dna)\b", header, re.I)
    if m:
        kind = m.group(1).lower()
        if kind == 'protein':
            return sbol3.SBO_PROTEIN
        if kind == 'rna':
            return sbol3.SBO_RNA
        return sbol3.SBO_DNA

    # Heuristic: if only AAs (20 std + *), treat as protein
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY\*\s]+", seq.upper()):
        return sbol3.SBO_PROTEIN

    # Simple heuristic: default DNA
    return sbol3.SBO_DNA

def guess_roles(header: str, filename: str) -> List[str]:
    """Return a list of SO role CURIEs; explicit role=SO:... wins."""
    # Allow explicit roles: role=SO:0000167 (multiple allowed)
    roles = re.findall(r"\brole=(SO:\d{7})\b", header, re.I)
    if roles:
        return [r.upper() for r in roles]

    hay = f"{header} | {os.path.basename(filename)}".lower()
    for pattern, so_curie in KEYWORD_TO_ROLE:
        if re.search(pattern, hay):
            return [so_curie]

    # Fallback: if looks like protein, assign no DNA roles
    return []

def clean_display_id(name: str) -> str:
    """Make a machine-safe displayId (no spaces, URL-safe-ish)."""
    # Keep letters, numbers, _ and -
    did = re.sub(r"[^A-Za-z0-9_\-\.]", "_", name.strip())
    # Avoid starting with non-letter
    if not re.match(r"^[A-Za-z_]", did):
        did = f"x_{did}"
    return did

def parse_header_tags(header: str) -> Dict[str, str]:
    """Extract simple key=value tags from header (type=, role=, host=, etc.)."""
    tags = dict(re.findall(r"(\b[a-zA-Z_]+)=([^\s|]+)", header))
    return tags

# ---------- FASTA parsing ----------

def read_fasta(path: str) -> List[Tuple[str, str]]:
    """Return list of (header, sequence) for all records in a FASTA file."""
    records = []
    header = None
    seq_lines: List[str] = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines).replace(" ", "").replace("\t", "")))
                header = line[1:].strip()
                seq_lines = []
            else:
                if line.strip():
                    seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines).replace(" ", "").replace("\t", "")))
    return records

def record_name_from_header(header: str) -> str:
    """Take the first token before space or pipe as human name."""
    # Prefer anything before first pipe, else before first space
    if "|" in header:
        return header.split("|", 1)[0].strip()
    return header.split(None, 1)[0].strip()

# ---------- SBOL document construction ----------

def make_component(doc: sbol3.Document,
                   name: str,
                   display_id: str,
                   seq_str: str,
                   comp_type: str,
                   roles: List[str],
                   extra_ann: Optional[Dict[str, str]] = None) -> sbol3.Component:
    """Create a Component with an attached Sequence and add to doc."""
    c = sbol3.Component(display_id, comp_type)
    c.name = name
    if roles:
        c.roles = roles

    # Sequence
    seq_id = f"{display_id}_seq"
    # Choose encoding based on type
    if comp_type == sbol3.SBO_PROTEIN:
        encoding = sbol3.IUPAC_PROTEIN_ENCODING
    else:
        encoding = sbol3.IUPAC_DNA_ENCODING

    s = sbol3.Sequence(seq_id, elements=seq_str, encoding=encoding)
    doc.add(s)
    c.sequences = [s]

    # Optional annotations (free-text via 'description'; more structured
    # annotations typically use custom namespaces; keep simple here)
    if extra_ann:
        desc_fields = [f"{k}={v}" for k, v in sorted(extra_ann.items()) if k not in {'type', 'role'}]
        if desc_fields:
            c.description = " | ".join(desc_fields)

    doc.add(c)
    return c

def build_document(inputs: List[str],
                   collection_name: str,
                   collection_id: Optional[str]) -> sbol3.Document:
    sbol3.set_namespace("https://example.org/parts/")  # customize if you have a real namespace
    doc = sbol3.Document()

    # Create a Collection
    col_id = collection_id or clean_display_id(collection_name or "parts_collection")
    coll = sbol3.Collection(col_id)
    coll.name = collection_name or "Parts Collection"
    doc.add(coll)

    # Iterate files & records
    for path in inputs:
        for header, seq in read_fasta(path):
            if not seq:
                continue
            name = record_name_from_header(header)
            display_id = clean_display_id(name)

            comp_type = guess_component_type(seq, header, path)
            roles = guess_roles(header, path)
            tags = parse_header_tags(header)

            c = make_component(
                doc=doc,
                name=name,
                display_id=display_id,
                seq_str=seq,
                comp_type=comp_type,
                roles=roles,
                extra_ann=tags
            )
            coll.members.append(c)

    return doc

# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA parts to SBOL3.")
    parser.add_argument("inputs", nargs="+", help="Input FASTA files (can be multiple).")
    parser.add_argument("--out", required=True, help="Output base name (no extension).")
    parser.add_argument("--collection-name", default="Parts Collection", help="SBOL Collection.name")
    parser.add_argument("--collection-id", default=None, help="SBOL Collection display_id (optional)")
    args = parser.parse_args()

    # Build doc
    doc = build_document(args.inputs, args.collection_name, args.collection_id)

    # Write both formats
    xml_path = f"{args.out}.xml"
    json_path = f"{args.out}.json"

    # Note: sbol3.Document.write infers format from extension
    doc.write(xml_path)
    doc.write(json_path)

    print(f"Wrote: {xml_path}")
    print(f"Wrote: {json_path}")

if __name__ == "__main__":
    main()