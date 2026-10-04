"""jplace version 3 serialisation."""

from __future__ import annotations

from .placement import EdgePlacement

FIELDS = ["edge_num", "likelihood", "like_weight_ratio",
          "distal_length", "pendant_length"]


def build_jplace(tree_with_edge_nums: str,
                 placements_by_query: dict[str, list[EdgePlacement]]) -> dict:
    """Assemble a jplace v3 document.

    'likelihood' is the optimised log-likelihood of the whole alignment;
    'like_weight_ratio' is the relative support of the edge (normalised over
    all edges of this query) -- not a probability of species correctness.
    """
    placements = []
    for qid, places in placements_by_query.items():
        placements.append({
            "p": [[p.edge_num, p.log_likelihood, p.weight,
                   p.distal_length, p.pendant_length] for p in places],
            "n": [qid],
        })
    return {
        "version": 3,
        "fields": list(FIELDS),
        "tree": tree_with_edge_nums,
        "placements": placements,
        "metadata": {
            "info": "eDNA placement on a fixed reference tree, JC69 model; "
                    "like_weight_ratio is relative support across edges, not "
                    "a probability of correct species assignment."
        },
    }
