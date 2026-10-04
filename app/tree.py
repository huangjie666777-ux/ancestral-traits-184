"""Unrooted binary tree graph with stable edge numbering and jplace Newick output."""

from __future__ import annotations

from dataclasses import dataclass


def quote_name(name: str) -> str:
    """Quote a Newick label so commas, parentheses, quotes and braces survive."""
    plain = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
                "0123456789_.-")
    if name and all(ch in plain for ch in name):
        return name
    return "'" + name.replace("'", "''") + "'"


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
            return quote_name(g.leaf_name[node])
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


@dataclass
class _NT:
    children: list
    name: str | None
    length: float | None
    edge_num: int | None


class _Tok:
    def __init__(self, text: str) -> None:
        self.s = text.strip()
        if self.s.endswith(";"):
            self.s = self.s[:-1]
        self.i = 0

    def skip_ws(self) -> None:
        while self.i < len(self.s) and self.s[self.i].isspace():
            self.i += 1

    def label(self) -> str:
        self.skip_ws()
        if self.i >= len(self.s) or self.s[self.i] != "'":
            start = self.i
            while self.i < len(self.s) and self.s[self.i] not in ",:;()[]{}":
                self.i += 1
            return self.s[start:self.i].strip()
        self.i += 1
        out: list[str] = []
        while True:
            if self.i >= len(self.s):
                raise ValueError("unterminated quoted label")
            ch = self.s[self.i]
            if ch == "'":
                if self.i + 1 < len(self.s) and self.s[self.i + 1] == "'":
                    out.append("'")
                    self.i += 2
                else:
                    self.i += 1
                    return "".join(out)
            else:
                out.append(ch)
                self.i += 1

    def expect(self, ch: str) -> None:
        self.skip_ws()
        if self.i >= len(self.s) or self.s[self.i] != ch:
            raise ValueError("expected '" + ch + "' at position " + str(self.i))
        self.i += 1

    def number(self) -> float:
        self.skip_ws()
        start = self.i
        while self.i < len(self.s) and (self.s[self.i].isdigit()
                                        or self.s[self.i] in ".eE+-"):
            self.i += 1
        token = self.s[start:self.i].strip()
        if not token:
            raise ValueError("missing number at position " + str(start))
        return float(token)

    def integer(self) -> int:
        self.skip_ws()
        start = self.i
        while self.i < len(self.s) and self.s[self.i].isdigit():
            self.i += 1
        token = self.s[start:self.i]
        if not token:
            raise ValueError("missing edge number at position " + str(start))
        return int(token)

    def branch(self) -> tuple[float | None, int | None]:
        length = None
        edge_num = None
        self.skip_ws()
        if self.i < len(self.s) and self.s[self.i] == ":":
            self.i += 1
            length = self.number()
        self.skip_ws()
        if self.i < len(self.s) and self.s[self.i] == "{":
            self.i += 1
            edge_num = self.integer()
            self.expect("}")
        return length, edge_num


def parse_annotated_newick(text: str) -> dict:
    """Parse jplace Newick (quoted labels, ':length{edge_num}') into geometry.

    Returns {'lengths': {edge_num: length},
             'subtree': {edge_num: (child_edges, leaf_name|None)},
             'edge_order': [edge_num, ...]}.
    Edges are oriented parent -> child; 'subtree' lists the immediate
    child edges on the distal side (not the edge itself, not deeper
    descendants), so distal subtree mass can be accumulated exactly once.
    """
    tok = _Tok(text)

    def node() -> _NT:
        tok.skip_ws()
        children: list[_NT] = []
        if tok.i < len(tok.s) and tok.s[tok.i] == "(":
            tok.i += 1
            while True:
                children.append(node())
                tok.skip_ws()
                if tok.i < len(tok.s) and tok.s[tok.i] == ",":
                    tok.i += 1
                    continue
                break
            tok.expect(")")
        name = tok.label() or None
        length, edge_num = tok.branch()
        return _NT(children, name, length, edge_num)

    root = node()
    tok.skip_ws()
    if tok.i != len(tok.s):
        raise ValueError("trailing characters at position " + str(tok.i))

    lengths: dict[int, float] = {}
    child_of: dict[int, list[int]] = {}
    leaf_of: dict[int, str | None] = {}
    edge_order: list[int] = []

    def walk(nd: _NT) -> None:
        for child in nd.children:
            walk(child)
        if nd.edge_num is not None:
            if nd.length is None:
                raise ValueError("edge " + str(nd.edge_num) + " has no length")
            if nd.edge_num in lengths:
                raise ValueError("duplicate edge number " + str(nd.edge_num))
            lengths[nd.edge_num] = float(nd.length)
            child_of[nd.edge_num] = [c.edge_num for c in nd.children
                                     if c.edge_num is not None]
            leaf_of[nd.edge_num] = nd.name if not nd.children else None
            edge_order.append(nd.edge_num)

    walk(root)
    return {"lengths": lengths, "subtree": child_of,
            "leaf": leaf_of, "edge_order": edge_order}
