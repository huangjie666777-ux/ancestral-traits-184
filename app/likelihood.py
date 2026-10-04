"""JC69 model and Felsenstein pruning with per-column scaling.

All conditional-likelihood messages are computed for every directed edge of
the fixed reference tree (two-pass rerooting), scaled per alignment column to
avoid underflow on long alignments.  Grafting a query onto an edge then only
needs the two messages of that edge plus the query leaf likelihoods.
"""

from __future__ import annotations

import math

import numpy as np

from .tree import TreeGraph

BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}
PI = 0.25


def jc69_matrix(t: float) -> np.ndarray:
    """Transition probability matrix P(t) under JC69."""
    a = math.exp(-4.0 * t / 3.0)
    p_same = 0.25 + 0.75 * a
    p_diff = 0.25 - 0.25 * a
    m = np.full((4, 4), p_diff)
    np.fill_diagonal(m, p_same)
    return m


def obs_vector(compatible: frozenset) -> np.ndarray:
    """Observation likelihood over states; gaps/unknown give no evidence."""
    v = np.zeros(4)
    if not compatible:  # gap: no base evidence
        v[:] = 1.0
    else:
        for b in compatible:
            v[BASE_INDEX[b]] = 1.0
    return v


class DirectedMessages:
    """Per-directed-edge conditional likelihood messages with log scales.

    msg[eid][endpoint] = (L, log_scale) where L has shape (n_cols, 4): the
    conditional likelihood of the data on the far side of the edge given the
    state at that endpoint, divided by per-column scaling factors.
    """

    def __init__(self, g: TreeGraph, leaf_obs: dict[int, np.ndarray]):
        self.msg: list[dict[int, tuple[np.ndarray, np.ndarray]]] = [
            {} for _ in range(len(g.edges))]
        self._compute(g, leaf_obs)

    def _compute(self, g: TreeGraph, leaf_obs: dict[int, np.ndarray]) -> None:
        n_cols = next(iter(leaf_obs.values())).shape[0]
        incoming: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}

        def combine(node: int, skip: int | None):
            # Message leaving 'node', combining all incoming except from 'skip'.
            if node in leaf_obs and len(g.adj[node]) == 1:
                return leaf_obs[node], np.zeros(n_cols)
            out = np.ones((n_cols, 4))
            scale = np.zeros(n_cols)
            for nbr, eid in g.adj[node]:
                if nbr == skip:
                    continue
                m, s = incoming[(node, nbr)]
                p = jc69_matrix(g.edges[eid].length)
                out *= m @ p.T
                scale += s
            mx = out.max(axis=1)
            mx[mx == 0.0] = 1.0
            return out / mx[:, None], scale + np.log(mx)

        # postorder from an arbitrary root, then preorder (rerooting)
        root = next(iter(g.adj))
        parent: dict[int, int] = {root: -1}
        order = [root]
        for u in order:
            for v, _ in g.adj[u]:
                if v not in parent:
                    parent[v] = u
                    order.append(v)
        for u in reversed(order):
            if u == root:
                continue
            incoming[(parent[u], u)] = combine(u, parent[u])
        for u in order[1:]:
            incoming[(u, parent[u])] = combine(parent[u], u)

        for eid, e in enumerate(g.edges):
            self.msg[eid][e.u] = incoming[(e.u, e.v)]
            self.msg[eid][e.v] = incoming[(e.v, e.u)]


def graft_column_loglik(msg_u, msg_v, len_u, len_v, query_obs, pendant):
    """Log-likelihood per column with the query grafted on the edge.

    len_u/len_v split the original edge length (both >= 0, sum preserved);
    pendant is the new leaf branch length in [0, 2].  Returns
    (loglik_per_column, total_loglik).
    """
    lu, su = msg_u
    lv, sv = msg_v
    # msg_u arrives at u from the v side, so it travels to the join point
    # along the v-side segment (len_v), and symmetrically for msg_v.
    a = lu @ jc69_matrix(len_v).T           # v-side message at the join point
    b = lv @ jc69_matrix(len_u).T           # u-side message at the join point
    q = query_obs @ jc69_matrix(pendant).T  # query leaf likelihood at join point
    col = PI * (a * b * q).sum(axis=1)
    col = np.maximum(col, 1e-300)
    logl = np.log(col) + su + sv
    return logl, float(logl.sum())
