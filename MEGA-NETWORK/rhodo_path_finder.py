from classes_functions import *
from kegg_parser import parse_kgml_to_tsv
import xml.etree.ElementTree as ET
from sys import argv
import os
import requests

PATHWAY_FOLDER = 'Pathways'
os.makedirs(PATHWAY_FOLDER, exist_ok=True)

print('\nHERO Bioinformatics Tool\nHigh Performance Rhodococcus Opacus\nby Team BACMAN\n')
print('You can use this tool to find the pathway that connects your initial compound to your desired product.\n')

print("\nThe pathways currently available in this local directory are:")
available_files = [f for f in os.listdir(PATHWAY_FOLDER) if f.endswith('.kgml')]
available_pathways = set(f.removesuffix('.kgml') for f in available_files)
for p in sorted(available_pathways):
    print(f" - {p.replace('_',' ')}")

response = input("\nWould you like to download more? (y/n) ").lower()
while response not in ['y','n']:
    response = input("Invalid response. Program only accepts Y or N.\nWould you like to download more? (y/n) ").lower()

if response == 'y':
    while True:
        pway_request = input("Type the pathway code (e.g. hsa00010) or press enter to finish: ").strip()
        if not pway_request:
            break
        if not (len(pway_request) == 8 and pway_request[0:3].islower() and pway_request[3:].isdigit()):
            print("Invalid KEGG pathway code. Try again.")
            continue
        try:
            print(f"Fetching {pway_request} from KEGG...")
            url = f"https://rest.kegg.jp/get/{pway_request}/kgml"
            r = requests.get(url, timeout=10)
            r.raise_for_status()
            try:
                root = ET.fromstring(r.text)
                pway_name = root.attrib.get("title", "unknown_pathway").replace(" ", "_") 
            except :
                pway_name = pway_request
            local_kgml = os.path.join(PATHWAY_FOLDER, f"{pway_name}.kgml")
            with open(local_kgml, 'w', encoding='utf-8') as f:
                f.write(r.text)
            print(f"Downloaded {pway_request} to {local_kgml}...\nParsing...")
            parse_kgml_to_tsv(local_kgml)
            print("Parsing completed.")
        except Exception as e:
            print(f"Failed to fetch or parse {pway_request}: {e}")


# Build all Graphs
edge_files = [f for f in os.listdir(PATHWAY_FOLDER) if f.endswith('_edges.tsv')]
for file in edge_files:
    MetabolicGraph.from_tsv(file)

# Build MEGAGRAPH
megagraph = MetabolicGraph.megaGraph()

node_list = list(Node.nodes.keys())
print('\n\nSELECTING THE STARTING AND FINAL COMPOUND')
while True:
        start = input('\nChoose a starting compound: ').strip()
        if start.lower() in node_list:
             break
        matches = [n for n in node_list if start.lower() in n]
        if matches:
            print(f"\nNo exact match for '{start}'. Did you mean:")
            for i, m in enumerate(matches[:5], 1):
                print(f"  {i}. {m}")
            print("Type one of those exactly, or press enter to input a different compound.")
        else:
            print(f"Compound '{start}' not found in the parsed .tsv files. You can select one from the available compounds found in the .tsv files in the Pathways folder.")

while True:
        end = input('\nChoose a final compound: ').strip()
        if end.lower() in node_list:
            break
        matches = [n for n in node_list if end.lower() in n]
        if matches:
            print(f"\nNo exact match for '{end}'. Did you mean:")
            for i, m in enumerate(matches[:5], 1):
                print(f"  {i}. {m}")
            print("Type one of those exactly, or press enter to input a different compound.")
        else:
            print(f"Compound '{end}' not found in the parsed .tsv files. You can select one from the available compounds found in the .tsv files in the Pathways folder.")


# Find shortest path
print('\n\nPATHWAY FINDING\n')
print(f"Looking for the shortest path from '{start}' to '{end}'...")
dist, path, ecs = megagraph.dijkstra(start, end)
if path is None:
    print('No path exists among the current nodes.')
    print('\nRetroPath 2.0 will be used to predict the pathway.')
else:
    print('\nPathway Steps:')
    for i in range(0,dist):
        print(f'{i+1}. {path[i].name} ---[{ecs[i]}]---> {path[i+1].name}')

# Check which enzymes you already have and which ones you need
# Build plasmid
print('\n\nPLASMID BUILDING\n')
print('Coming soon...')

print('\nThank you for using HERO :)\n')