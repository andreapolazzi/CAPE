from pathway_functions import view, download
from classes_functions import MetabolicGraph
from retropath_converter import build_retopath_network
# from retropath2_wrapper import RetroPath2 # only works in Linux
import sys, os, difflib
from colorama import init, Fore, Style
init()

PATHWAY_FOLDER = 'Pathways'
os.makedirs(PATHWAY_FOLDER, exist_ok=True)
args = sys.argv[1:]

print(Fore.GREEN +
      "\nHERO Bioinformatics Tool\n"
      "H" + Style.RESET_ALL + "igh p" +
      Fore.GREEN + "E" + Style.RESET_ALL +
      "rformance " + Fore.GREEN +
      "R" + Style.RESET_ALL +
      "hodococcus " + Fore.GREEN +
      "O" + Style.RESET_ALL +
      "pacus\nby Team BACMAN\n")
print("You can use this tool to find the pathway that connects your initial "
      "compound to your desired product\nand to build the plasmid that you can "
      "insert into your Rhodococcus (coming soon).")

if len(args) == 1 and args[0] == 'pathways':
    view()

elif len(args) == 1 and args[0] == 'download':
    download()

elif len(args) == 2:
    # 1) Build mega graph
    if not os.listdir(PATHWAY_FOLDER):
        print('Warning: the Pathways folder is missing or empty. '
              'Download pathways before attempting to find shortest paths.')
        sys.exit(1)

    for fname in os.listdir(PATHWAY_FOLDER):
        if not fname.endswith('_edges.tsv'):
            continue
        MetabolicGraph.from_tsv(fname)

    MetabolicGraph.mega = MetabolicGraph.megaGraph()

    # 2) Select start and end
    node_list = list(n.lower() for n in MetabolicGraph.mega.nodes)
    start, end = args
    start = start.strip()
    end = end.strip()

    def prompt_choice(name, node_list):
        while name.lower() not in node_list:
            # fuzzy matches
            matches = [n for n in node_list if name.lower() in n]
            matches += difflib.get_close_matches(name.lower(), node_list, n=10, cutoff=0.1)
            matches = sorted(set(matches))
            if matches:
                print(f"\nNo exact match for '{name}' in the current nodes. Did you mean:")
                for i, m in enumerate(matches, 1):
                    print(f"  {i}. {m}")
                name = input("Type one of those exactly, or type 'RP' to have RetroPath 2.0 predict the path, or press Enter to retry: ").strip()
            else:
                name = input(f"'{name}' not found. Try again: ").strip()
        return name

    start = prompt_choice(start, node_list)
    end   = prompt_choice(end,   node_list)

    # 3) Connect start and/or end to central metabolites if they are not in the nodes
    if start == 'RP':
        # Connenct start to central metabolites through RetroPath
        pass
    if end == 'RP':
        # Connenct end to central metabolites through RetroPath
        pass

    # 4) Find up to 3 shortest paths
    print(Fore.GREEN + '\nPATHWAY FINDING\n' + Style.RESET_ALL)
    print(f"Looking for the top 3 shortest paths from '{start}' to '{end}'...")
    all_paths = MetabolicGraph.mega.shortest_paths(start, end)

    # 5) Fallback to RetroPath if path can't be found (not yet fully implemented)
    if not all_paths:
        print('No path exists among the current nodes.')
        print('\nRetroPath 2.0 will be used to fill in the gaps and predict the pathway.')

        # Connect start to metaboliti principali
        # create source file using the start compound
        # create sink file using the metaboliti principali
        # create rules file using the inchis of the compounds
        # r_code = RetroPath2(sink_file='/path/to/sink/file',source_file='/path/to/source/file',
            # rules_file='/path/to/rules/file',outdir='start_RP_results.csv')
        
        # Connect end to metaboliti principali
        # create source file using the end compound
        # create sink file using the metaboliti principali
        # create rules file using the inchis of the compounds
        # r_code = RetroPath2(sink_file='/path/to/sink/file',source_file='/path/to/source/file',
            # rules_file='/path/to/rules/file',outdir='end_RP_results.csv')

        # Build graphs from RP results and join them to mega graph
        build_retopath_network('start_RP_results.csv')
        build_retopath_network('end_RP_results.csv')
        start_pway_graph = MetabolicGraph.from_tsv('start_RP_results_edges.tsv')
        end_pway_graph = MetabolicGraph.from_tsv('end_RP_results_edges.tsv')
        start_pway_graph.join2megaGraph()
        end_pway_graph.join2megaGraph()

        # Retry finding the shortest path
        all_paths = MetabolicGraph.mega.shortest_paths(start, end)
        if not all_paths:
            print("RetroPath 2.0 couldn't find a path.")
            sys.exit(1)

    # 6) Print each candidate path
    WEIGHT_COLOR = {1: Fore.GREEN, 10: Fore.YELLOW, 100: Fore.RED}
    for idx, (dist, path_nodes, path_ecs, path_rxns) in enumerate(sorted(all_paths, key=lambda x: x[0]), start=1):
        print(f"\n{Fore.CYAN}Path #{idx}{Style.RESET_ALL}")
        for i in range(len(path_nodes) - 1):
            u, v = path_nodes[i], path_nodes[i+1]
            ec = path_ecs[i]
            # find the actual edge object
            chosen = next(
                (e for e in u.outgoing
                if v in e.products and (e.ec == ec or e.name == ec)),
                None
            )
            if not chosen:
                chosen = next((e for e in u.outgoing if v in e.products), None)
            w = chosen.weight if chosen else 1
            col = WEIGHT_COLOR.get(w, Fore.WHITE)
            print(f"  {i+1}. {u.name} "
                f"{Fore.WHITE}---[{Style.RESET_ALL}"
                f"{col}{ec.replace(' ', ' or ')}{Style.RESET_ALL}"
                f"{Fore.WHITE}]--->{Style.RESET_ALL} "
                f"{v.name}")

    # 7) Prompt user to choose one
    while True:
        choice = input(f"\nSelect which path to use for plasmid building [1-{len(all_paths)}]: ")
        if choice.isdigit() and 1 <= (n := int(choice)) <= len(all_paths):
            break
        print("Invalid selection, please choose a number in range.")

    dist, path_nodes, path_ecs, path_rxns = all_paths[n-1]

    # 8) Visualize only the chosen path
    sub = MetabolicGraph.mega.get_subgraph(path_nodes, path_ecs, path_rxns)
    sub.visualize(f"Path #{n} from {start} to {end}", start, end)

    # 9) Placeholder for plasmid building
    print(Fore.GREEN + '\n\nPLASMID BUILDING\n' + Style.RESET_ALL)
    print('Coming soon...')
    print('\n\nThank you for using HERO :)\n')

else:
    print("\nUsage:\n"
          "py rhodo_path_finder.py pathways            # List available pathways\n"
          "py rhodo_path_finder.py download            # Fetch new pathways\n"
          "py rhodo_path_finder.py <start> <end>       # Find and visualize shortest path from start to end compound")
