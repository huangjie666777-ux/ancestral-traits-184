"""Per-edge placement optimisation and like-weight ratios."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .likelihood import DirectedMessages, graft_column_loglik
from .tree import TreeGraph

MAX_PENDANT = 2.0


@dataclass
class EdgePlacement:
    edge_num: int
    log_likelihood: float
    distal_length: float    # join point -> child-side node (Edge.v)
    pendant_length: float   # new leaf branch length
    weight: float = 0.0     # like_weight_ratio, normalised over all edges


def _objective(dm: DirectedMessages, g: TreeGraph, eid: int,
               query_obs: np.ndarray, len_u: float, pendant: float) -> float:
    e = g.edges[eid]
    len_v = e.length - len_u
    if len_u < -1e-12 or len_v < -1e-12 or not (0.0 <= pendant <= MAX_PENDANT):
        return -math.inf
    _, total = graft_column_loglik(dm.msg[eid][e.u], dm.msg[eid][e.v],
                                   max(len_u, 0.0), max(len_v, 0.0),
                                   query_obs, pendant)
    return total


def optimise_edge(dm: DirectedMessages, g: TreeGraph, eid: int,
                  query_obs: np.ndarray) -> EdgePlacement:
    """Maximise the grafted-tree log-likelihood over split point and pendant
    length for one edge (coarse grid + pattern-search refinement)."""
    e = g.edges[eid]
    best = (-math.inf, 0.0, 0.0)
    for len_u in np.linspace(0.0, e.length, 11):
        for pendant in np.linspace(0.0, MAX_PENDANT, 11):
            ll = _objective(dm, g, eid, query_obs, len_u, pendant)
            if ll > best[0]:
                best = (ll, float(len_u), float(pendant))
    ll, len_u, pendant = best
    step_u = max(e.length / 10.0, 1e-4)
    step_t = MAX_PENDANT / 10.0
    while max(step_u, step_t) > 1e-6:
        improved = False
        for du, dt in ((step_u, 0), (-step_u, 0), (0, step_t), (0, -step_t),
                       (step_u, step_t), (-step_u, -step_t)):
            nu = min(max(len_u + du, 0.0), e.length)
            nt = min(max(pendant + dt, 0.0), MAX_PENDANT)
            cand = _objective(dm, g, eid, query_obs, nu, nt)
            if cand > ll + 1e-12:
                ll, len_u, pendant = cand, nu, nt
                improved = True
        if not improved:
            step_u *= 0.5
            step_t *= 0.5
    return EdgePlacement(edge_num=eid, log_likelihood=ll,
                         distal_length=e.length - len_u,
                         pendant_length=pendant)


def place_query(dm: DirectedMessages, g: TreeGraph,
                query_obs: np.ndarray) -> list[EdgePlacement]:
    """Optimise every edge, normalise weights, sort by likelihood then edge."""
    placements = [optimise_edge(dm, g, eid, query_obs)
                  for eid in range(len(g.edges))]
    top = max(p.log_likelihood for p in placements)
    raw = [math.exp(p.log_likelihood - top) for p in placements]
    total = sum(raw)
    for p, r in zip(placements, raw):
        p.weight = r / total
    placements.sort(key=lambda p: (-p.log_likelihood, p.edge_num))
    return placements


def has_base_evidence(states: list[frozenset]) -> bool:
    """True if at least one column gives discriminating base evidence."""
    return any(0 < len(s) < 4 for s in states)
