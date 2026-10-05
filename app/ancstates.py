"""Ancestral state reconstruction by weighted Sankoff (Fitch) parsimony.

The reference geometry is read from a stored placement job's jplace tree
(edge-numbered Newick): leaves and edge numbers are preserved exactly, and
the unrooted tree is re-rooted at the internal node adjacent to the outgroup
leaf.  Each edge costs one transition lookup parent-state -> child-state in a
direction-aware integer cost matrix (None forbids that direction); branch
lengths are not used and no intermediate states are chained along an edge.

Global (whole-tree) optima only: down/up (reroot) messages decide which
states are possible at each node and which ordered state pairs are possible
on each edge *among globally optimal histories* -- never a single arbitrary
optimum and never the Cartesian product of the two endpoint state sets.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .tree import parse_annotated_newick

INF = float("inf")


class AncestralError(ValueError):
    """A located problem that rejects the whole reconstruction request."""


@dataclass
class RTNode:
    key: object                 # str leaf id, or frozenset of incident edge nums
    leaf: str | None = None
    children: list["RTNode"] = field(default_factory=list)
    parent: "RTNode | None" = None
    edge_to_parent: int | None = None   # edge_num on the parent link

    @property
    def is_leaf(self) -> bool:
        return self.leaf is not None


def build_rooted_tree(tree_text: str, outgroup: str):
    """Build an undirected graph from jplace Newick and root it at the
    internal node adjacent to the outgroup leaf.

    Returns (nodes_by_key, root, edge_endpoints, leaf_ids) where
    edge_endpoints maps edge_num -> (parent-side key, child-side key) after
    rooting.  All leaves and original edge numbers are retained; the outgroup
    pendant edge is the root -> outgroup edge.
    """
    geo = parse_annotated_newick(tree_text)
    child_of = geo["subtree"]
    leaf_of = geo["leaf"]
    all_edges = set(geo["lengths"])

    appears_as_child: set[int] = set()
    for kids in child_of.values():
        appears_as_child.update(kids)
    root_edges = frozenset(e for e in all_edges if e not in appears_as_child)

    def parent_group(e: int) -> frozenset:
        if e in root_edges:
            return root_edges
        for p, kids in child_of.items():
            if e in kids:
                return frozenset([p] + kids)
        raise AncestralError("edge " + str(e) + " has no parent-side node")

    def child_key(e: int):
        if leaf_of[e] is not None:
            return leaf_of[e]
        return frozenset([e] + child_of[e])

    adj: dict[object, list[tuple[object, int]]] = {}
    leaf_edges: dict[str, int] = {}
    for e in all_edges:
        a, b = parent_group(e), child_key(e)
        adj.setdefault(a, []).append((b, e))
        adj.setdefault(b, []).append((a, e))
        if leaf_of[e] is not None:
            leaf_edges[leaf_of[e]] = e

    if outgroup not in leaf_edges:
        raise AncestralError("outgroup leaf '" + outgroup + "' is not a leaf "
                             "of the reference tree")
    root_key = parent_group(leaf_edges[outgroup])
    if not isinstance(root_key, frozenset):
        raise AncestralError("node adjacent to outgroup '" + outgroup +
                             "' is not an internal node")

    nodes: dict[object, RTNode] = {key: RTNode(key) for key in adj}
    for key, node in nodes.items():
        if isinstance(key, str):
            node.leaf = key

    order = [root_key]
    parent_of: dict[object, object] = {root_key: None}
    edge_between: dict[object, int] = {}
    for key in order:
        for nb, e in adj[key]:
            if nb in parent_of:
                continue
            parent_of[nb] = key
            edge_between[nb] = e
            order.append(nb)
    if len(order) != len(nodes):
        raise AncestralError("reference tree geometry is not connected")

    root = nodes[root_key]
    for key in order[1:]:
        node, pnode = nodes[key], nodes[parent_of[key]]
        node.parent = pnode
        node.edge_to_parent = edge_between[key]
        pnode.children.append(node)

    endpoints: dict[int, tuple[object, object]] = {}
    for key in order[1:]:
        endpoints[edge_between[key]] = (parent_of[key], key)
    return nodes, root, endpoints, sorted(leaf_edges)


def _feasible_sets(root: RTNode, allowed: dict[str, frozenset],
                   n_states: int, cost: list[list]) -> dict:
    """Bottom-up set of states from which each subtree is realizable."""
    feas: dict[object, frozenset] = {}
    order: list[RTNode] = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(u.children)
    for u in reversed(order):
        if u.is_leaf:
            feas[u.key] = allowed[u.leaf]
            continue
        acc = set(range(n_states))
        for v in u.children:
            fv = feas[v.key]
            can = {a for a in range(n_states)
                   if any(cost[a][b] is not None for b in fv)}
            acc &= can
        feas[u.key] = frozenset(acc)
    return feas


def reconstruct(tree_text: str, outgroup: str, states: list[str],
                allowed: dict[str, frozenset],
                cost: list[list]) -> dict:
    """Run the full reconstruction and return a JSON-serialisable report."""
    nodes, root, endpoints, leaf_ids = build_rooted_tree(tree_text, outgroup)
    k = len(states)

    # ---- feasibility diagnosis (no half results when infeasible) ----
    feas = _feasible_sets(root, allowed, k, cost)
    if not feas[root.key]:
        reasons = []
        for e in sorted(endpoints):
            _, ck = endpoints[e]
            distal = feas[ck]
            if not any(cost[a][b] is not None
                       for a in range(k) for b in distal):
                target = ("leaf:" + nodes[ck].leaf) if nodes[ck].is_leaf else                     "internal node " + _node_label(nodes[ck])
                reasons.append("edge " + str(e) + " towards " + target +
                               ": no permitted parent->child state pair can "
                               "reach the states feasible for that subtree")
        if not reasons:
            reasons.append("no joint node assignment satisfies the leaf "
                           "constraints under the permitted transitions")
        return {"feasible": False, "reasons": reasons}

    # ---- postorder Sankoff DP with optimal-history counts ----
    down: dict[object, list] = {}
    cnt: dict[object, list] = {}
    order: list[RTNode] = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(u.children)
    for u in reversed(order):
        if u.is_leaf:
            down[u.key] = [0 if i in allowed[u.leaf] else INF
                           for i in range(k)]
            cnt[u.key] = [1 if i in allowed[u.leaf] else 0 for i in range(k)]
            continue
        du, cu = [], []
        for s in range(k):
            total = 0.0
            ways = 1
            ok = True
            for v in u.children:
                best = INF
                for t in range(k):
                    if cost[s][t] is None:
                        continue
                    val = down[v.key][t] + cost[s][t]
                    if val < best:
                        best = val
                if best == INF:
                    ok = False
                    break
                total += best
                n = sum(cnt[v.key][t] for t in range(k)
                        if cost[s][t] is not None
                        and down[v.key][t] + cost[s][t] == best)
                ways *= n
            du.append(total if ok else INF)
            cu.append(ways if ok else 0)
        down[u.key], cnt[u.key] = du, cu

    opt = min(down[root.key])
    opt_count = sum(cnt[root.key][s] for s in range(k)
                    if down[root.key][s] == opt)

    # ---- reroot (up) messages: best cost outside each subtree ----
    up: dict[object, list] = {root.key: [0.0] * k}
    for u in order:
        for v in u.children:
            dv = down[v.key]
            best_in = [min((dv[t] + cost[s][t] for t in range(k)
                            if cost[s][t] is not None), default=INF)
                       for s in range(k)]
            sib = [down[u.key][s] - best_in[s] for s in range(k)]
            row = []
            for t in range(k):
                row.append(min((up[u.key][s] + sib[s] + cost[s][t]
                                for s in range(k) if cost[s][t] is not None),
                               default=INF))
            up[v.key] = row

    possible: dict[object, list] = {}
    for key in nodes:
        possible[key] = [s for s in range(k)
                         if down[key][s] + up[key][s] == opt]

    # ---- per-edge ordered parent->child pairs that occur in some optimum ----
    edge_reports = []
    for e in sorted(endpoints):
        pk, ck = endpoints[e]
        dv = down[ck]
        best_in = [min((dv[t] + cost[s][t] for t in range(k)
                        if cost[s][t] is not None), default=INF)
                   for s in range(k)]
        sib = [down[pk][s] - best_in[s] for s in range(k)]
        pairs = []
        any_same = any_diff = False
        for s in range(k):
            for t in range(k):
                if cost[s][t] is None:
                    continue
                total = up[pk][s] + sib[s] + cost[s][t] + dv[t]
                if total == opt:
                    pairs.append([states[s], states[t]])
                    if s == t:
                        any_same = True
                    else:
                        any_diff = True
        if any_diff and any_same:
            change = "possibly_change"
        elif any_diff:
            change = "always_change"
        else:
            change = "always_same"
        edge_reports.append({
            "edge_num": e,
            "parent_node": _node_label(nodes[pk]),
            "child_node": _node_label(nodes[ck]),
            "possible_pairs": pairs,
            "change": change,
        })

    node_reports = []
    for u in order:
        node_reports.append({
            "node_id": _node_label(u),
            "kind": "leaf" if u.is_leaf else "internal",
            "leaf_id": u.leaf,
            "edges": sorted(u.key) if isinstance(u.key, frozenset) else
                     [u.edge_to_parent],
            "is_root": u is root,
            "parent_node": _node_label(u.parent) if u.parent else None,
            "possible_states": [states[s] for s in possible[u.key]],
        })

    return {
        "feasible": True,
        "minimum_cost": opt,
        "optimal_histories": str(opt_count),
        "root_node": _node_label(root),
        "nodes": node_reports,
        "edges": edge_reports,
    }


def _node_label(u: RTNode) -> str:
    if u.leaf is not None:
        return "leaf:" + u.leaf
    return "internal:[" + ",".join(str(e) for e in sorted(u.key)) + "]"
