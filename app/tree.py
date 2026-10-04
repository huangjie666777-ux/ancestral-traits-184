"""Unrooted binary tree graph with stable edge numbering and jplace Newick output."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Edge:
    u: int          # parent-side node (orientation used for distal_length)
    v: int          # child-side node
    length: float


class TreeGraph:
    def __init__(self) -> None:
        self.adj: dict[int, list[tuple[int, int]]] = {}  # node -> [(neighbor, edge_id)]
        self.edges: list[Edge] = []
        self.leaf_name: dict[int, str] = {}

    def add_node(self) -> int:
        nid = len(self.adj)
        self.adj[nid] = []
        return nid

    def add_edge(self, u: int, v: int, length: float) -> int:
        eid = len(self.edges)
        self.edges.append(Edge(u, v, length))
        self.adj[u].append((v, eid))
        self.adj[v].append((u, eid))
        return eid

    @property
    def n_nodes(self) -> int:
        return len(self.adj)


def from_phylotree(tree) -> TreeGraph:
    """Convert a Bio.Phylo tree to an unrooted graph.

    A bifurcating root is dissolved: the two root edges merge into one edge
    whose length is the sum, keeping the unrooted topology and total lengths.
    """
    g = TreeGraph()
    node_of: dict[int, int] = {}

    def node_for(clade) -> int:
        key = id(clade)
        if key not in node_of:
            nid = g.add_node()
            node_of[key] = nid
            if clade.is_terminal():
                g.leaf_name[nid] = clade.name
        return node_of[key]

    def build(clade) -> int:
        nid = node_for(clade)
        for child in clade.clades:
            cid = build(child)
            g.add_edge(nid, cid, float(child.branch_length))
        return nid

    root = tree.root
    if len(root.clades) == 2:
        c1, c2 = root.clades
        n1, n2 = build(c1), build(c2)
        g.add_edge(n1, n2, float(c1.branch_length) + float(c2.branch_length))
    else:
        build(root)
    return g


def renumber_edges(g: TreeGraph) -> list[int]:
    """Assign stable edge numbers by DFS from the alphabetically first leaf.

    Returns a list mapping edge_id -> edge_num; edge objects are reordered so
    that edge_id == edge_num afterwards.
    """
    first_leaf = min(g.leaf_name, key=lambda n: g.leaf_name[n])
    order: list[int] = []
    visited = {first_leaf}

    def dfs(u: int) -> None:
        for v, eid in sorted(g.adj[u], key=lambda t: g.edges[t[1]].length):
            if eid in order:
                continue
            order.append(eid)
            if v not in visited:
                visited.add(v)
                dfs(v)

    dfs(first_leaf)
    remap = {eid: num for num, eid in enumerate(order)}
    g.edges = [g.edges[eid] for eid in order]
    g.adj = {n: [] for n in g.adj}
    for eid, e in enumerate(g.edges):
        g.adj[e.u].append((e.v, eid))
        g.adj[e.v].append((e.u, eid))
    return remap


def orient_and_serialize(g: TreeGraph) -> str:
    """Serialize with jplace edge-number annotations, rooted at a degree-3 node.

    Edges are oriented parent->child from that root; Edge.u/Edge.v are updated
    in place so 'distal_length' (toward v) is well defined.
    """
    root = next((n for n in g.adj if len(g.adj[n]) == 3),
                next(iter(g.adj)))
    visited = {root}

    def fmt(node: int) -> str:
        if node in g.leaf_name:
            return g.leaf_name[node]
        parts = []
        for v, eid in list(g.adj[node]):
            if v in visited:
                continue
            visited.add(v)
            e = g.edges[eid]
            e.u, e.v = node, v
            parts.append(fmt(v) + ":" + repr(float(e.length)) + "{" + str(eid) + "}")
        return "(" + ",".join(parts) + ")"

    body = fmt(root)
    return body + ";"
