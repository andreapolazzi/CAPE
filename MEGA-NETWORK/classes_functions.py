from __future__ import annotations
from itertools import count
import csv, heapq

class Node:
    def __init__(self, name: str, inchi: str, ingoing: list[Edge], outgoing: list[Edge], g: MetabolicGraph):
        self.name = name
        self.inchi = inchi
        self.ingoing = ingoing
        self.outgoing = outgoing
        self.org = [g.org] if g.org else []
        self.pathway = [g.pathway] if g.pathway else []

dont_traverse_list = ['h2o','h+','nad+','nadh','nadph','nadp+','co2','oxygen','adp','atp']

class Edge:
    def __init__(self, name: str, ec: str, source: Node, target: Node, g: MetabolicGraph, weight: float = 1):
        self.name = name
        self.ec = ec
        self.source = source
        self.target = target
        self.org = [g.org] if g.org else []
        self.pathway = [g.pathway] if g.pathway else []
        if source.name.lower() in dont_traverse_list or target.name.lower() in dont_traverse_list:
            self.weight = float('inf') # so that the path cannot traverse water but has to traverse a more important compound from the reaction step
        else:
            self.weight = weight
        key = f"{g.org}:{name}:{source.name}->{target.name}"
        source.outgoing.append(self)
        target.ingoing.append(self)

class MetabolicGraph:
    graphs: list[MetabolicGraph] = []
    mega: MetabolicGraph = None
    def __init__(self):
        self.nodes: dict[str, Node] = {}  # keyed by lowercase
        self.edges: dict[str, Edge] = {}
        self.org: str = ''
        self.pathway: str = ''
        MetabolicGraph.graphs.append(self)

    def add_node(self, name: str, inchi: str = ''):
        key = name.lower()
        if key not in self.nodes:
            # create Node with original name
            n = Node(name=name, inchi=inchi, ingoing=[], outgoing=[], g=self)
            self.nodes[key] = n
        else:
            n = self.nodes[key]
            if self.org and self.org not in n.org:
                n.org.append(self.org)
            if self.pathway and self.pathway not in n.pathway:
                n.pathway.append(self.pathway)


    def add_edge(self, enzyme: str, ec: str, source_name: str, target_name: str, source_inchi: str = '', target_inchi: str = '', weight: float = 1):
        # ensure nodes exist
        self.add_node(source_name, source_inchi)
        self.add_node(target_name, target_inchi)
        src = self.nodes[source_name.lower()]
        tgt = self.nodes[target_name.lower()]
        key = f"{self.org}:{enzyme}:{src.name}->{tgt.name}"
        if key in self.edges:
            e = self.edges[key]
            if self.org and self.org not in e.org:
                e.org.append(self.org)
            if self.pathway and self.pathway not in e.pathway:
                e.pathway.append(self.pathway)

        else:
            e = Edge(name=enzyme, ec=ec, source=src, target=tgt, g=self, weight=weight)
            self.edges[key] = e

    @classmethod
    def from_tsv(cls, tsv_file: str, weight: float = 1) -> MetabolicGraph:
        """Load a metabolic network from TSV with columns: enzyme_name, EC, source_node_name, target_node_name"""
        graph = cls()
        graph.org = tsv_file[0:3]
        graph.pathway = tsv_file[3:]
        with open(f"Pathways/{tsv_file}", 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                enz = row['enzyme_name']
                ec = row.get('EC', '')
                src = row['source']
                tgt = row['target']
                clr = row['color']
                if clr == 'green':
                    w = 1
                elif clr == 'white':
                    w = 10
                else: # retropath
                    w = 100
                graph.add_edge(enzyme=enz, ec=ec, source_name=src, target_name=tgt, weight=w )
        return graph

    @classmethod
    def megaGraph(cls) -> MetabolicGraph:
        """Combine all MetabolicGraph objects into one"""
        mega = cls()
        mega.org = 'MegaGraph'
        for g in cls.graphs:
            # add all nodes (org and pathway lists will update automatically, without duplicate nodes)
            for n in g.nodes.values():
                mega.add_node(n.name, n.inchi)
        for g in cls.graphs:
            # add all edges (org and pathway lists will update automatically, without duplicate edges)
            for e in g.edges.values():
                mega.add_edge(enzyme=e.name, ec=e.ec, source_name=e.source.name,
                              target_name=e.target.name, source_inchi=e.source.inchi,
                              target_inchi=e.target.inchi, weight=e.weight)
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
            mega.add_edge(enzyme=e.name, ec=e.ec,
                source_name=e.source.name, target_name=e.target.name,
                source_inchi=e.source.inchi, target_inchi=e.target.inchi,
                weight=e.weight)
        MetabolicGraph.mega = mega
        return mega

    def dijkstra(self, start: str, goal: str) -> tuple[float, list[Node], list[str]]:
        """Returns (distance, path_nodes, path_ec_numbers) or (None, None, None) if unreachable."""
        try:
            start_node = self.nodes[start.lower()]
            goal_node  = self.nodes[goal.lower()]
        except KeyError:
            print('nodes not found')
            return None, None, None

        dist: dict[Node, float] = {n: float('inf') for n in self.nodes.values()}
        prev: dict[Node, tuple[Node, Edge] | None] = {n: None for n in self.nodes.values()}
        dist[start_node] = 0

        heap: list[tuple[float, int, Node]] = []
        counter = count()
        heapq.heappush(heap, (0, next(counter), start_node))
        visited: set[Node] = set()

        while heap:
            curr_dist, _, u = heapq.heappop(heap)
            if u in visited:
                continue
            visited.add(u)

            if u is goal_node:
                break

            for e in u.outgoing:
                v = e.target
                if v in visited:
                    continue
                new_dist = curr_dist + e.weight
                if new_dist < dist[v]:
                    dist[v] = new_dist
                    prev[v] = (u, e)
                    heapq.heappush(heap, (new_dist, next(counter), v))

        if dist[goal_node] == float('inf'):
            print('infinite distance')
            return None, None, None

        # Reconstruct node‐path and collect EC numbers
        path_nodes: list[Node] = []
        path_ecs: list[str] = []
        node = goal_node
        while node:
            path_nodes.append(node)
            p = prev[node]
            if p is None:
                break
            prev_node, edge_used = p
            path_ecs.append(edge_used.ec or edge_used.name or "")
            node = prev_node
        path_nodes.reverse()
        path_ecs.reverse()

        return dist[goal_node], path_nodes, path_ecs