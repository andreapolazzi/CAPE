from classes_functions import *
from kegg_parser import parse_kgml_to_tsv
import xml.etree.ElementTree as ET
import sys
import os
import glob
import requests
from colorama import init, Fore, Back, Style
init()

sys.setrecursionlimit(10_000)

PATHWAY_FOLDER = 'Pathways'
os.makedirs(PATHWAY_FOLDER, exist_ok=True)
GRAPH_VAR_PATH = os.path.join(PATHWAY_FOLDER, 'graph_variable')
graph_variable_present = os.path.exists(GRAPH_VAR_PATH)

print(Fore.GREEN + '\nHERO Bioinformatics Tool\nH' + Style.RESET_ALL + 'igh p' + Fore.GREEN + 'E' + Style.RESET_ALL + 'rformance ' + Fore.GREEN + 'R' + Style.RESET_ALL + 'hodococcus ' + Fore.GREEN + 'O' + Style.RESET_ALL + 'pacus\nby Team BACMAN\n')
print('You can use this tool to find the pathway that connects your initial compound to your desired product\nand to build the plasmid that you can insert into your Rhodococcus (coming soon).\n')

if graph_variable_present:
    MetabolicGraph.mega = MetabolicGraph.load_graph_variable()
else:
    edge_files = [f for f in os.listdir(PATHWAY_FOLDER) if f.endswith('_edges.tsv')]
    for file in edge_files:
        MetabolicGraph.from_tsv(file)
    MetabolicGraph.megaGraph().save_graph_variable()

print("\nThe pathways currently available in the mega-network are:")
available_files = [f for f in os.listdir(PATHWAY_FOLDER) if f.endswith('.kgml')]
available_pathways = set(f.removesuffix('.kgml') for f in available_files)
for p in sorted(available_pathways):
    print(f" - {p.replace('_',' ')}")

response = input("\nWould you like to download more? (y/n) ").lower()
while response not in ['y','n']:
    response = input("Invalid response. Program only accepts Y or N.\nWould you like to download more? (y/n) ").lower()

if response.lower() == 'y':
    while True:
        line = input("Type one or multiple space-separated pathway codes (e.g. hsa00010 hsa00020),\nor type the organism code only (e.g. hsa) to download all its pathways,\nor press Enter to quit downloading: ").strip()
        if not line:
            break
        elif len(line) == 3:
            response = requests.get(f"https://rest.kegg.jp/link/pathway/{line}")
            if response.status_code == 200:
                pathway_ids = set() # set to avoid duplicates
                for line in response.text.strip().split("\n"):
                    pathway_id = line.split("\t")[1].removeprefix("path:")
                    pathway_ids.add(pathway_id)
                line = ' '.join(sorted(list(pathway_ids)))
            else:
                print(f"Organism code '{line}' not found by KEGG API.")

        codes = line.split()
        for pway_request in codes:
            if not (len(pway_request) == 8 and pway_request[:3].islower() and pway_request[3:].isdigit()):
                print(f"Invalid KEGG pathway code '{pway_request}'. Skipping.")
                continue

            try:
                kgml_filename = f"{pway_request}.kgml"
                if glob.glob(f"Pathways/*{pway_request}*"):
                    print(f"\nPathway {pway_request} already downloaded. Skipping.")
                    continue

                print(f"\nFetching {pway_request} from KEGG…")
                url = f"https://rest.kegg.jp/get/{pway_request}/kgml"
                r = requests.get(url, timeout=10)
                r.raise_for_status()

                try:
                    root = ET.fromstring(r.text)
                    raw_title = root.attrib.get("title", pway_request)
                    kgml_title = raw_title.replace("/", "and").replace(" ", "_")
                except Exception:
                    kgml_title = pway_request

                kgml_path = os.path.join(PATHWAY_FOLDER, f"{kgml_title}.kgml")
                with open(kgml_path, 'w', encoding='utf-8') as f:
                    f.write(r.text)
                print(f"Saved as {kgml_title}.kgml.\nParsing…")
                parse_kgml_to_tsv(kgml_path)
                print("Parsing completed.")

                new_pathway_graph = MetabolicGraph.from_tsv(f'{pway_request}_edges.tsv')
                new_pathway_graph.join2megaGraph()
                print("Joined parsed pathway to mega-network.")

            except Exception as e:
                print(f"Failed to fetch or parse '{pway_request}': {e}\n")
                # Clean files if they were created
                if kgml_path and os.path.exists(kgml_path):
                    os.remove(kgml_path)
                if f'{pway_request}_nodes.tsv' and f'{pway_request}_edges.tsv' and os.path.exists(f'{pway_request}_nodes.tsv') and os.path.exists(f'{pway_request}_edges.tsv'):
                    os.remove(f'{pway_request}_nodes.tsv')
                    os.remove(f'{pway_request}_edhges.tsv')

    MetabolicGraph.mega.save_graph_variable() # save new mega graph after new pathways were added

node_list = list(MetabolicGraph.mega.nodes.keys())
print(Fore.GREEN + '\n\nSELECTING THE STARTING AND FINAL COMPOUND' + Style.RESET_ALL)
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
print(Fore.GREEN + '\n\nPATHWAY FINDING\n' + Style.RESET_ALL)
print(f"Looking for the shortest path from '{start}' to '{end}'...")
dist, path, ecs = MetabolicGraph.mega.dijkstra(start, end)
if path is None:
    print('No path exists among the current nodes.')
    print('\nRetroPath 2.0 will be used to predict the pathway.')
else:
    print('\nPathway Steps:')
    for i in range(dist):
        arrow_tip = Fore.GREEN + ']--->' + Style.RESET_ALL
        arrow_end = Fore.GREEN + '---[' + Style.RESET_ALL
        print(f"{i+1}. {path[i].name} {arrow_end}{ecs[i].replace(' ', ' or ')}{arrow_tip} {path[i+1].name}")

# Check which enzymes you already have and which ones you need
# Build plasmid
print(Fore.GREEN + '\n\nPLASMID BUILDING\n' + Style.RESET_ALL)
print('Coming soon...')

print('\n\nThank you for using HERO :)\n')