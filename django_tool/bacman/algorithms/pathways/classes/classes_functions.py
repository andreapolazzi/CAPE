from __future__ import annotations
from itertools import count
import csv, heapq
import matplotlib.pyplot as plt
import pydot
from PIL import Image
from io import BytesIO
from copy import deepcopy


# List of compounds to avoid in traversal
dont_traverse_list = [
    # Energy carriers
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump",
    "ctp", "cdp", "cmp", "ppi", "pi", "inorganic phosphate", "orthophosphate", "diphosphate", "reduced fmn",
    # Electron carriers
    "nad+", "nadh", "nadp+", "nadph", "fad", "fadh2", "fmn", "fmnh2","electron-transferring flavoprotein",
    "reduced acceptor", "acceptor", "[reduced nadph---hemoprotein reductase]", "[oxidized nadph---hemoprotein reductase]",
    "reduced ferredoxin", "oxidized ferredoxin", "ferredoxin",
    # Coenzymes and group carriers
    "coenzyme a", "coa", "acetyl-coa", "s-adenosylmethionine", "sam", "thf", "tetrahydrofolate",
    "biotin", "lipoic acid", "pyridoxal phosphate", "plp", "tpp", "thiamine pyrophosphate", "acyl-carrier protein",
    # Gases and small byproducts
    "co2", "h2o", "nh3", "h+", "o2", "oxygen", "h2o2", "ch4", "formate", "formaldehyde",
    # Sugars and nucleotides in transfer reactions
    "udp-glucose", "cdp-choline", "cmp-neu5ac", "prpp", "ribose-5-phosphate", "dump", "dtmp",
    # Miscellaneous
    "glutathione", "gsh", "gssg", "carnitine", "creatine", "phosphocreatine",
    "ubiquinone", "ubiquinol", "cytochrome c", "ferricytochrome c", "fe2+", "fe3+",
]




class Node:
    def __init__(self, name: str, inchi: str, g: MetabolicGraph):
        self.name = name
        self.inchi = inchi
        self.ingoing: list[Edge] = []
        self.outgoing: list[Edge] = []
        self.org: list[str] = [g.org] if g.org else []
        self.pathway: list[str] = [g.pathway] if g.pathway else []




class Edge:
    def __init__(self, name: str, ec: str, reagents: list[Node], products: list[Node], reaction_id: str, g: MetabolicGraph, weight: float = 1):
        self.name = name
        self.ec = ec
        self.reagents = reagents
        self.products = products
        self.reaction_id = reaction_id
        self.org: list[str] = [g.org] if g.org else []
        self.pathway: list[str] = [g.pathway] if g.pathway else []
        # determine weight; penalize edges with undesired compounds
        self.weight = weight
        for r in self.reagents + self.products:
            if r.name.lower() in dont_traverse_list:
                self.weight = float('inf')
        # register this edge with nodes
        for r in reagents:
            r.outgoing.append(self)
        for p in products:
            p.ingoing.append(self)




class MetabolicGraph:
    graphs: list[MetabolicGraph] = []
    mega: MetabolicGraph = None


    def __init__(self):
        self.nodes: dict[str, Node] = {}
        self.edges: dict[str, Edge] = {}
        self.org: str = ''
        self.pathway: str = ''
        MetabolicGraph.graphs.append(self)


    def add_node(self, name: str, inchi: str = '') -> Node:
        key = name.lower()
        if key not in self.nodes:
            node = Node(name=name, inchi=inchi, g=self)
            self.nodes[key] = node
        else:
            node = self.nodes[key]
            if self.org and self.org not in node.org:
                node.org.append(self.org)
            if self.pathway and self.pathway not in node.pathway:
                node.pathway.append(self.pathway)
        return node


    def add_edge(self,
                 enzyme: str, ec: str,
                 reaction_id: str,
                 reagent_names: list[str], product_names: list[str],
                 reagent_inchis: list[str] = None, product_inchis: list[str] = None,
                 weight: float = 1):
        reagent_inchis = reagent_inchis or [''] * len(reagent_names)
        product_inchis = product_inchis or [''] * len(product_names)

        reagents = [self.add_node(n, i) for n, i in zip(reagent_names, reagent_inchis)]
        products = [self.add_node(n, i) for n, i in zip(product_names, product_inchis)]

        # unique reaction key
        key = f"{self.org}:{enzyme}:{'+'.join(r.name for r in reagents)}->" \
              f"{'+'.join(p.name for p in products)}"
        if key in self.edges:
            edge = self.edges[key]
            if self.org and self.org not in edge.org:
                edge.org.append(self.org)
            if self.pathway and self.pathway not in edge.pathway:
                edge.pathway.append(self.pathway)
        else:
            edge = Edge(name=enzyme, ec=ec, reagents=reagents, products=products, reaction_id=reaction_id, g=self, weight=weight)
            self.edges[key] = edge


    @classmethod
    def from_tsv(cls, tsv_file: str) -> MetabolicGraph:
        graph = cls()
        graph.org = tsv_file[:3]
        graph.pathway = tsv_file[3:]
        with open(f"/work/igem/iGEM/django_tool/bacman/algorithms/pathways/Pathways/{tsv_file}", 'r', newline='', encoding='utf-8') as f:
        #with open(f"./bacman/algorithms/pathways/Pathways/{tsv_file}", 'r', newline='', encoding='utf-8') as f: # to run in local
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                enz = row['enzyme_name']
                ec = row.get('EC', '')
                reaction_id = row['reaction']
                src = row['source']
                tgt = row['target']
                clr = row['color']
                weight = 1 if clr=='green' else 10 if clr=='white' else 100
                graph.add_edge(
                    enzyme=enz,
                    ec=ec,
                    reaction_id=reaction_id,
                    reagent_names=[src],
                    product_names=[tgt],
                    reagent_inchis=[''],
                    product_inchis=[''],
                    weight=weight
                )
        return graph


    @classmethod
    def megaGraph(cls) -> MetabolicGraph:
        mega = cls()
        mega.org = 'MegaGraph'
        for g in cls.graphs:
            for n in g.nodes.values():
                mega.add_node(n.name, n.inchi)
        for g in cls.graphs:
            for e in g.edges.values():
                mega.add_edge(
                    enzyme=e.name,
                    ec=e.ec,
                    reaction_id=e.reaction_id,
                    reagent_names=[r.name for r in e.reagents],
                    product_names=[p.name for p in e.products],
                    reagent_inchis=[r.inchi for r in e.reagents],
                    product_inchis=[p.inchi for p in e.products],
                    weight=e.weight
                )
        cls.mega = mega
        return mega


    def join2megaGraph(self) -> MetabolicGraph:
        """Join a graph to the mega graph"""
        if MetabolicGraph.mega is None:
            MetabolicGraph.mega = MetabolicGraph()
            MetabolicGraph.mega.org = 'MegaGraph'
        mega = MetabolicGraph.mega
        for n in self.nodes.values():
            mega.add_node(n.name, n.inchi)
        for e in self.edges.values():
            mega.add_edge(
                    enzyme=e.name,
                    ec=e.ec,
                    reaction_id=e.reaction_id,
                    reagent_names=[r.name for r in e.reagents],
                    product_names=[p.name for p in e.products],
                    reagent_inchis=[r.inchi for r in e.reagents],
                    product_inchis=[p.inchi for p in e.products],
                    weight=e.weight
                )
        MetabolicGraph.mega = mega
        return mega


    # Dijkstra (one path)
    def dijkstra(self, start: str, goal: str, banned_nodes: set = None, banned_edges: set = None) -> tuple[float, list[Node], list[str], list[str]]:
        """Single‐path Dijkstra from `start` to `goal`, skipping any node in banned_nodes or any directed edge (u,v) in banned_edges.
        Returns (distance, [nodes], [edge_ECs]) or (None, None, None) if no path."""
        banned_nodes = banned_nodes or set()
        banned_edges = banned_edges or set()

        try:
            start_node = self.nodes[start.lower()]
            goal_node  = self.nodes[goal.lower()]
        except KeyError:
            print('nodes not found')
            return None, None, None, None

        dist = {n: float('inf') for n in self.nodes.values()}
        prev = {n: None for n in self.nodes.values()}
        dist[start_node] = 0.0

        heap = []
        tie = count()
        heapq.heappush(heap, (0.0, next(tie), start_node))
        visited = set()

        while heap:
            curr_dist, _, u = heapq.heappop(heap)
            if u in visited:
                continue
            visited.add(u)
            if u in banned_nodes:
                continue
            if u is goal_node:
                break

            for e in u.outgoing:
                if e in banned_edges:
                    continue
                for v in e.products:
                    if v in banned_nodes or v in visited:
                        continue
                    new_dist = curr_dist + e.weight
                    if new_dist < dist[v]:
                        dist[v] = new_dist
                        prev[v] = (u, e)
                        heapq.heappush(heap, (new_dist, next(tie), v))

        if dist[goal_node] == float('inf'):
            return None, None, None, None

        # Reconstruct path nodes, ECs, and rxns
        path_nodes, path_ecs, path_rxns = [], [], []
        node = goal_node
        while node:
            path_nodes.append(node)
            p = prev[node]
            if p is None:
                break
            u, e = p
            path_ecs.append(e.ec or e.name)
            path_rxns.append(e.reaction_id)
            node = u

        path_nodes.reverse()
        path_ecs.reverse()
        path_rxns.reverse()
        return dist[goal_node], path_nodes, path_ecs, path_rxns


    # Yen’s K-shortest simple paths using Dijkstra
    def shortest_paths(self, start: str, goal: str, K: int = 3) -> list[tuple[float, list[Node], list[str]]]:
        """Yen’s algorithm for up to K simple shortest paths from start to goal.
        Returns a list of (distance, [nodes], [ECs]) tuples."""
        # Get first shortest path
        first = self.dijkstra(start, goal, set(), set())
        if first[0] is None:
            return []
        paths = [first]  # accepted paths
        candidates = []  # min‑heap of (dist, tie, nodes, ecs, rxns, banned_nodes, banned_edges)
        tie = count()

        # Build up to K paths
        for k in range(1, K):
            _, prev_nodes, prev_ecs, prev_rxns = paths[k-1]
            for i in range(len(prev_nodes)-1):
                root_nodes = prev_nodes[:i+1]
                root_ecs   = prev_ecs[:i]
                root_rxns  = prev_rxns[:i]
                spur_node  = root_nodes[-1]

                # Ban exactly those Edge objects used in earlier A‑paths with same root
                banned_edges = set()
                for (_, nodesA, ecsA, rxnsA) in paths:
                    if nodesA[:i+1] == root_nodes:
                        u = nodesA[i]
                        v = nodesA[i+1]
                        rxn = rxnsA[i]

                        # Find the unique Edge in u.outgoing matching that step
                        for e in u.outgoing:
                            if rxn == e.reaction_id and v in e.products:
                                banned_edges.add(e)
                                break

                # Ban the root nodes except spur
                banned_nodes = set(root_nodes[:-1])

                spur = self.dijkstra(spur_node.name,
                                    goal,
                                    banned_nodes=banned_nodes,
                                    banned_edges=banned_edges)
                if spur[0] is None:
                    continue

                spur_dist, spur_nodes, spur_ecs, spur_rxns = spur

                # Compute root distance
                root_dist = 0.0
                for j, (u, v) in enumerate(zip(root_nodes, root_nodes[1:])):
                    desired_rxn = root_rxns[j]
                    for e in u.outgoing:
                        if v in e.products and e.reaction_id == desired_rxn:
                            root_dist += e.weight
                            break

                total_dist  = root_dist + spur_dist
                total_nodes = root_nodes[:-1] + spur_nodes
                total_ecs   = root_ecs   + spur_ecs
                total_rxns  = root_rxns  + spur_rxns

                heapq.heappush(candidates, (
                    total_dist,
                    next(tie),
                    total_nodes,
                    total_ecs,
                    total_rxns,
                    frozenset(banned_nodes),
                    frozenset(banned_edges)
                ))

            # Select the next best path
            while candidates:
                dist_c, _, nodes_c, ecs_c, rxns_c, _, _ = heapq.heappop(candidates)
                if not any(nodes_c == p[1] for p in paths):
                    paths.append((dist_c, nodes_c, ecs_c, rxns_c))
                    break
            else:
                break

        return paths


    def get_subgraph(self,path: list[Node],ecs:  list[str],rxns: list[str]) -> MetabolicGraph:
        """Build a subgraph containing every full reaction (all reagents & products)
        for each reaction_id in 'rxns'."""
        subgraph = MetabolicGraph()
        subgraph.org     = self.org
        subgraph.pathway = f"{path[0].name}_to_{path[-1].name}"
        for rxn_id in rxns:
            # Get all edges in the mega‑graph that correspond to that reaction
            for e in self.edges.values():
                if e.reaction_id != rxn_id:
                    continue
                # Now add every reagent→product pair for this reaction
                for r in e.reagents:
                    for p in e.products:
                        subgraph.add_edge(
                            enzyme         = e.name,
                            ec             = e.ec,
                            reaction_id    = e.reaction_id,
                            reagent_names  = [r.name],
                            product_names  = [p.name],
                            reagent_inchis = [r.inchi],
                            product_inchis = [p.inchi],
                            weight         = e.weight
                        )
        return subgraph


    @staticmethod
    def directional_add_edge(graph, src, dst, label):
        edge_fwd = graph.get_edge(src, dst)
        edge_rev = graph.get_edge(dst, src)
        # Case: already reversible -> do nothing
        if any(e.get("dir") == "both" for e in edge_fwd + edge_rev):
            return
        # Case: both directions exist -> make one bidirectional edge
        if edge_fwd and edge_rev:
            for e in edge_fwd + edge_rev:
                graph.del_edge(e.get_source(), e.get_destination())
            all_labels = set(e.get("label") for e in edge_fwd + edge_rev)
            all_labels.add(label)
            graph.add_edge(pydot.Edge(src, dst, color="black", dir="both",
                                    label=", ".join(sorted(all_labels)), fontsize="8"))
        # Case: reverse exists -> add forward
        elif edge_rev:
            graph.add_edge(pydot.Edge(src, dst, color="black", label=label, fontsize="8"))
        # Case: forward exists -> add reverse
        elif edge_fwd:
            graph.add_edge(pydot.Edge(dst, src, color="black", label=label, fontsize="8"))
        # Case: no edge yet -> add this one
        else:
            graph.add_edge(pydot.Edge(src, dst, color="black", label=label, fontsize="8"))


    def visualize(self, title: str, start: str, end: str):
        graph = pydot.Dot(graph_type='digraph', rankdir='LR')
        nodes_added = set()

        for e in self.edges.values():
            enzyme_label = f"{e.name} [EC:{e.ec}] (RN{e.reaction_id})"
            node_id = f"enz_{e.reaction_id}"

            # Reagent nodes
            for r in e.reagents:
                if r.name not in nodes_added:
                    fill = '#A2C8E2' if r.name.lower() in (start, end) else '#C8A2C8'
                    node = pydot.Node(r.name,shape='box',style='filled,rounded',fillcolor=fill)
                    graph.add_node(node)
                    nodes_added.add(r.name)

            # Enzyme nodes
            # Compute effective weights (ignoring infinite weights)
            weights = [e2.weight for e2 in self.edges.values() if e2.reaction_id == e.reaction_id and e2.weight != float('inf')]
            w = min(weights) if weights else float('inf')

            # Now map weight to a color
            if w == 1:
                fill = '#90ee90'
            elif w == 10:
                fill = '#ffff66'
            elif w == 100:
                fill = '#ff6666'
            else:
                fill = '#cccccc'

            enzyme_node = pydot.Node(node_id,label=enzyme_label,shape='box',style='filled,rounded',fillcolor=fill)
            graph.add_node(enzyme_node)
            nodes_added.add(node_id)

            # Product Nodes
            for p in e.products:
                if p.name not in nodes_added:
                    fill = '#A2C8E2' if p.name.lower() in (start, end) else '#C8A2C8'
                    node = pydot.Node(p.name,
                                    shape='box',
                                    style='filled,rounded',
                                    fillcolor=fill)
                    graph.add_node(node)
                    nodes_added.add(p.name)

            # Draw edges
            for r in e.reagents:
                MetabolicGraph.directional_add_edge(graph,src=r.name,dst=node_id,label='')
            for p in e.products:
                MetabolicGraph.directional_add_edge(graph,src=node_id,dst=p.name,label='')

        # Show image with plt
        png_str = graph.create_png()
        img = Image.open(BytesIO(png_str))
        plt.figure(figsize=(14, 10))
        plt.imshow(img)
        plt.axis('off')
        plt.title(title)
        plt.show()


