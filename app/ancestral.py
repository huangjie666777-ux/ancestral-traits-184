"""Ancestral state reconstruction on a stored placement job's reference tree.

The tree is re-rooted at the internal node adjacent to a user-chosen outgroup
leaf; every leaf and every original jplace edge number is preserved.  Each
directed edge (parent state -> child state) pays exactly one transition cost
from the user-supplied matrix (never multiplied by branch length, and no
intermediate states are chained along a single edge).  A Sankoff dynamic
program finds the minimum total cost over all node assignments consistent
with the per-leaf allowed-state sets (a set expresses uncertainty, not
simultaneous polymorphism; an all-states set marks an unknown trait).

The report covers *every* globally optimal history: per-node states that
occur in at least one optimum, per-edge ordered parent->child state pairs
that occur in at least one optimum (not the cartesian product of the two
endpoint state sets), and the exact number of optimal full assignments as a
decimal string.  Ambiguity is reported as possibility, never as probability.
"""

from __future__ import annotations

from .tree import parse_annotated_newick


class AncestralError(ValueError):
    """A located problem that rejects the whole reconstruction request."""


def _check_states(raw) -> list:
    if not isinstance(raw, list):
        raise AncestralError("'states' must be a list of 2..5 unique names")
    if not 2 <= len(raw) <= 5:
        raise AncestralError("expected 2..5 states, got " + str(len(raw)))
    seen: set = set()
    for i, name in enumerate(raw):
        if not isinstance(name, str) or not name:
            raise AncestralError("states[" + str(i) + "] must be a non-empty "
                                 "string")
        if name in seen:
            raise AncestralError("duplicate state name '" + name +
                                 "' at states[" + str(i) + "]")
        seen.add(name)
    return list(raw)


def _check_matrix(raw, states: list) -> list:
    n = len(states)
    if not isinstance(raw, list) or len(raw) != n:
        raise AncestralError("'transition_costs' must be a " + str(n) + "x" +
                             str(n) + " matrix in the order of 'states'")
    matrix: list = []
    for i, row in enumerate(raw):
        if not isinstance(row, list) or len(row) != n:
            raise AncestralError("transition_costs[" + str(i) + "] must be a "
                                 "list of " + str(n) + " entries")
        out_row: list = []
        for j, cell in enumerate(row):
            where = ("transition_costs[" + str(i) + "][" + str(j) + "] (" +
                     states[i] + " -> " + states[j] + ")")
            if isinstance(cell, bool):
                raise AncestralError(where + ": booleans are not allowed; use "
                                     "a non-negative integer or null")
            if cell is None:
                if i == j:
                    raise AncestralError(where + ": the diagonal must be 0, "
                                         "not null")
                out_row.append(None)
                continue
            if not isinstance(cell, int) or cell < 0:
                raise AncestralError(where + ": expected a non-negative "
                                     "integer or null, got " + repr(cell))
            if i == j and cell != 0:
                raise AncestralError(where + ": the diagonal must be 0, got "
                                     + str(cell))
            out_row.append(cell)
        matrix.append(out_row)
    return matrix


def _check_leaf_sets(raw, states: list, leaf_ids: set) -> dict:
    if not isinstance(raw, dict):
        raise AncestralError("'leaf_states' must be an object mapping every "
                             "reference leaf id to its allowed-state set")
    state_set = set(states)
    for key in raw:
        if key not in leaf_ids:
            raise AncestralError("unknown leaf '" + str(key) + "': not a leaf "
                                 "of the reference tree")
    missing = sorted(leaf_ids - set(raw))
    if missing:
        raise AncestralError("leaf '" + missing[0] + "' has no allowed-state "
                             "set; leaves must correspond exactly to the "
                             "reference tree (use the full state list to mark "
                             "an unknown trait)")
    allowed: dict = {}
    for leaf, values in raw.items():
        if not isinstance(values, list) or not values:
            raise AncestralError("leaf '" + leaf + "': allowed states must be "
                                 "a non-empty list")
        seen: set = set()
        for value in values:
            if not isinstance(value, str) or value not in state_set:
                raise AncestralError("leaf '" + leaf + "': unknown state '" +
                                     str(value) + "'")
            if value in seen:
                raise AncestralError("leaf '" + leaf + "': duplicate state '" +
                                     value + "' in its allowed set")
            seen.add(value)
        allowed[leaf] = seen
    return allowed


def _build_nodes(geo: dict) -> list:
    """Node list from parsed jplace geometry.

    Each edge's distal node is one node; the serialization root is one more.
    Returns [{'edges': sorted incident edge nums, 'leaf': name|None}].
    """
    subtree = geo["subtree"]
    referenced = set()
    for kids in subtree.values():
        referenced.update(kids)
    top = sorted(set(geo["lengths"]) - referenced)
    nodes = [{"edges": sorted([e] + subtree[e]), "leaf": geo["leaf"][e]}
             for e in geo["edge_order"]]
    nodes.append({"edges": top, "leaf": None})
    return nodes


def _reroot(nodes: list, outgroup: str):
    """Orient every edge parent->child from the node next to the outgroup leaf.

    Returns (root_idx, children, parent_edge, edge_endpoints) where
    edge_endpoints[e] = (parent_node_idx, child_node_idx).
    """
    og_leaf = next(i for i, nd in enumerate(nodes)
                   if nd["leaf"] == outgroup)
    og_edge = nodes[og_leaf]["edges"][0]
    root = next(i for i, nd in enumerate(nodes)
                if i != og_leaf and og_edge in nd["edges"])

    edge_to_nodes: dict = {}
    for i, nd in enumerate(nodes):
        for e in nd["edges"]:
            edge_to_nodes.setdefault(e, []).append(i)

    children: dict = {i: [] for i in range(len(nodes))}
    parent_edge: dict = {i: None for i in range(len(nodes))}
    edge_endpoints: dict = {}
    stack = [root]
    visited = {root}
    while stack:
        u = stack.pop()
        for e in nodes[u]["edges"]:
            v = next(w for w in edge_to_nodes[e] if w != u)
            if v in visited:
                continue
            visited.add(v)
            children[u].append((e, v))
            parent_edge[v] = e
            edge_endpoints[e] = (u, v)
            stack.append(v)
    return root, children, parent_edge, edge_endpoints


def _solve(nodes, root, children, allowed, matrix, states):
    """Sankoff DP.  Returns (min_cost, total_count, cost, count) or None."""
    k = len(states)
    cost: dict = {}
    count: dict = {}
    order: list = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(v for _, v in children[u])
    for u in reversed(order):
        leaf = nodes[u]["leaf"]
        cu: list = []
        cn: list = []
        for s in range(k):
            if leaf is not None:
                ok = states[s] in allowed[leaf]
                cu.append(0 if ok else None)
                cn.append(1 if ok else 0)
                continue
            total = 0
            ways = 1
            feasible = True
            for _, v in children[u]:
                best = None
                best_ways = 0
                for t in range(k):
                    w = matrix[s][t]
                    if w is None or cost[v][t] is None:
                        continue
                    val = cost[v][t] + w
                    if best is None or val < best:
                        best, best_ways = val, count[v][t]
                    elif val == best:
                        best_ways += count[v][t]
                if best is None:
                    feasible = False
                    break
                total += best
                ways *= best_ways
            cu.append(total if feasible else None)
            cn.append(ways if feasible else 0)
        cost[u] = cu
        count[u] = cn
    finite = [c for c in cost[root] if c is not None]
    if not finite:
        return None
    min_cost = min(finite)
    total_count = sum(count[root][s] for s in range(k)
                      if cost[root][s] == min_cost)
    return min_cost, total_count, cost, count


def _possible(nodes, root, children, cost, matrix, min_cost, k):
    """Downward pass: states/pairs realised by at least one global optimum."""
    possible: dict = {root: {s for s in range(k)
                             if cost[root][s] == min_cost}}
    edge_pairs: dict = {}
    stack = [root]
    while stack:
        u = stack.pop()
        for e, v in children[u]:
            poss_v: set = set()
            pairs: set = set()
            for s in possible[u]:
                best = None
                for t in range(k):
                    w = matrix[s][t]
                    if w is None or cost[v][t] is None:
                        continue
                    val = cost[v][t] + w
                    if best is None or val < best:
                        best = val
                for t in range(k):
                    w = matrix[s][t]
                    if w is not None and cost[v][t] is not None \
                            and cost[v][t] + w == best:
                        poss_v.add(t)
                        pairs.add((s, t))
            possible[v] = poss_v
            edge_pairs[e] = pairs
            stack.append(v)
    return possible, edge_pairs


def _node_id(node: dict) -> str:
    if node["leaf"] is not None:
        return "leaf:" + node["leaf"]
    return "internal:" + ",".join(str(e) for e in node["edges"])


def reconstruct(payload, jobs: dict) -> dict:
    """Validate the request and reconstruct ancestral states read-only."""
    if not isinstance(payload, dict):
        raise AncestralError("request body must be an object")
    job_id = payload.get("job_id")
    if not isinstance(job_id, str) or not job_id:
        raise AncestralError("'job_id' must be a non-empty string")
    if job_id not in jobs:
        raise AncestralError("unknown job_id '" + job_id + "'")
    states = _check_states(payload.get("states"))

    geo = parse_annotated_newick(jobs[job_id]["jplace"]["tree"])
    nodes = _build_nodes(geo)
    leaf_ids = {nd["leaf"] for nd in nodes if nd["leaf"] is not None}
    allowed = _check_leaf_sets(payload.get("leaf_states"), states, leaf_ids)

    outgroup = payload.get("outgroup")
    if not isinstance(outgroup, str) or not outgroup:
        raise AncestralError("'outgroup' must be a reference leaf id")
    if outgroup not in leaf_ids:
        raise AncestralError("unknown outgroup leaf '" + outgroup +
                             "': not a leaf of the reference tree")
    matrix = _check_matrix(payload.get("transition_costs"), states)

    root, children, parent_edge, edge_endpoints = _reroot(nodes, outgroup)
    solved = _solve(nodes, root, children, allowed, matrix, states)
    if solved is None:
        return {
            "job_id": job_id,
            "states": states,
            "feasible": False,
            "reason": ("no feasible history: the leaf constraints and the "
                       "forbidden (null) transitions admit no complete node "
                       "assignment; at least one required parent->child "
                       "transition is forbidden"),
        }
    min_cost, total_count, cost, count = solved
    possible, edge_pairs = _possible(nodes, root, children, cost, matrix,
                                     min_cost, len(states))

    node_reports = []
    stack = [root]
    while stack:
        u = stack.pop()
        nd = nodes[u]
        entry = {
            "node_id": _node_id(nd),
            "kind": "leaf" if nd["leaf"] is not None else "internal",
            "adjacent_edges": nd["edges"],
            "is_root": u == root,
            "parent_edge": parent_edge[u],
            "parent": (_node_id(nodes[edge_endpoints[parent_edge[u]][0]])
                       if parent_edge[u] is not None else None),
            "possible_states": sorted((states[s] for s in possible[u]),
                                      key=states.index),
        }
        if nd["leaf"] is not None:
            entry["leaf_id"] = nd["leaf"]
        node_reports.append(entry)
        stack.extend(v for _, v in children[u])

    edge_reports = []
    for e in sorted(edge_endpoints):
        u, v = edge_endpoints[e]
        pairs = sorted(edge_pairs[e])
        changed = [a != b for a, b in pairs]
        if all(changed):
            change = "always_change"
        elif any(changed):
            change = "may_change"
        else:
            change = "never_change"
        edge_reports.append({
            "edge_num": e,
            "parent": _node_id(nodes[u]),
            "child": _node_id(nodes[v]),
            "possible_pairs": [[states[a], states[b]] for a, b in pairs],
            "change": change,
        })

    return {
        "job_id": job_id,
        "states": states,
        "feasible": True,
        "root": _node_id(nodes[root]),
        "min_cost": min_cost,
        "optimal_histories": str(total_count),
        "nodes": node_reports,
        "edges": edge_reports,
    }
