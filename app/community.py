"""Community composition analysis: placement mass and tree KR (p=1) distances.

Each sample maps query ids to non-negative read counts.  A query's reads are
spread over *all* candidate placement locations using like_weight_ratio (not
only the best edge, and not collapsed onto the distal tip).  The resulting
measure lives on the fixed reference tree edges; the Kantorovich-Rubinstein
distance (Wasserstein-1, p=1) of two samples is the integral over every edge
of the absolute cumulative-mass difference in the distal subtree along that
edge, summed over edges.
"""

from __future__ import annotations

from .tree import parse_annotated_newick


class CommunityError(ValueError):
    """A located problem that rejects the whole community request."""


def _validate_samples(raw_samples) -> None:
    if not isinstance(raw_samples, list):
        raise CommunityError("'samples' must be a list")
    if not 2 <= len(raw_samples) <= 10:
        raise CommunityError("expected 2..10 samples, got " + str(len(raw_samples)))
    seen: set[str] = set()
    for idx, sample in enumerate(raw_samples):
        where = "sample[" + str(idx) + "]"
        if not isinstance(sample, dict):
            raise CommunityError(where + " must be an object")
        sid = sample.get("sample_id")
        if not isinstance(sid, str) or not sid:
            raise CommunityError(where + ".sample_id must be a non-empty string")
        if sid in seen:
            raise CommunityError("duplicate sample_id '" + sid + "'")
        seen.add(sid)
        job_id = sample.get("job_id")
        if not isinstance(job_id, str) or not job_id:
            raise CommunityError(where + " ('" + sid + "').job_id must be a "
                                 "non-empty string")
        counts = sample.get("counts")
        if not isinstance(counts, dict):
            raise CommunityError(where + " ('" + sid + "').counts must be an "
                                 "object mapping query id to read count")
        for qid, count in counts.items():
            if not isinstance(qid, str) or not qid:
                raise CommunityError(where + " ('" + sid + "'): query id keys "
                                     "must be non-empty strings")
            ok = isinstance(count, int) and not isinstance(count, bool) and count >= 0
            if not ok:
                raise CommunityError(where + " ('" + sid + "') query '" + qid +
                                     "': count must be a non-negative integer, got "
                                     + repr(count))


def validate_request(payload) -> list[dict]:
    if not isinstance(payload, dict) or "samples" not in payload:
        raise CommunityError("request body must contain 'samples'")
    _validate_samples(payload["samples"])
    return payload["samples"]


def resolve_jobs(samples: list[dict], jobs: dict) -> list[dict]:
    """Look up every referenced job and demand identical annotated Newick."""
    resolved: list[dict] = []
    reference_tree: str | None = None
    for sample in samples:
        job_id = sample["job_id"]
        if job_id not in jobs:
            raise CommunityError("sample '" + sample["sample_id"] + "': unknown "
                                 "job_id '" + job_id + "'")
        tree_text = jobs[job_id]["jplace"]["tree"]
        if reference_tree is None:
            reference_tree = tree_text
        elif tree_text != reference_tree:
            raise CommunityError(
                "sample '" + sample["sample_id"] + "' job '" + job_id + "' uses a "
                "different annotated reference tree (edge-numbered Newick must be "
                "identical to the first job's tree)")
        resolved.append(jobs[job_id])
    return resolved


def placements_by_id(job: dict) -> dict[str, list]:
    return {p["n"][0]: [dict(zip(job["jplace"]["fields"], row))
                        for row in p["p"]]
            for p in job["jplace"]["placements"]}


def build_measures(samples: list[dict], jobs: list[dict]):
    """Construct per-sample edge mass events from counts and placement weights.

    Returns (geometry, measures, report) where measures[i] maps
    edge_num -> sorted [(position_from_parent, normalized_mass), ...].
    Mass at a placement location is reads * like_weight_ratio for that query,
    accumulated over every query; zero counts are simply absent.
    """
    geo = parse_annotated_newick(jobs[0]["jplace"]["tree"])
    lengths = geo["lengths"]
    edge_set = set(lengths)

    measures: list[dict[int, list[tuple[float, float]]]] = []
    report: list[dict] = []

    for sample, job in zip(samples, jobs):
        places = placements_by_id(job)
        unplaceable = job.get("unplaceable", {})
        raw_events: dict[int, list[tuple[float, float]]] = {
            eid: [] for eid in edge_set}
        excluded: list[dict] = []
        valid_total = 0
        counts = sample["counts"]

        for qid in sorted(counts):
            count = counts[qid]
            if qid not in places and qid not in unplaceable:
                raise CommunityError("sample '" + sample["sample_id"] + "': "
                                     "unknown query id '" + qid + "'")
            if count == 0:
                continue
            if qid in places:
                candidates = places[qid]
                weight_sum = sum(c["like_weight_ratio"] for c in candidates)
                if not candidates or abs(weight_sum - 1.0) > 1e-7:
                    excluded.append({"id": qid, "count": count,
                                     "reason": "placement weights are not "
                                               "normalised in stored job"})
                    continue
                bad = [c for c in candidates
                       if c["edge_num"] not in edge_set
                       or not (0.0 <= c["distal_length"]
                               <= lengths[c["edge_num"]] + 1e-9)]
                if bad:
                    excluded.append({"id": qid, "count": count,
                                     "reason": "placement references an edge or "
                                               "position outside the reference "
                                               "tree"})
                    continue
                valid_total += count
                for c in candidates:
                    eid = c["edge_num"]
                    length = lengths[eid]
                    pos = min(max(length - c["distal_length"], 0.0), length)
                    raw_events[eid].append((pos, count * c["like_weight_ratio"]))
            else:
                excluded.append({"id": qid, "count": count,
                                 "reason": unplaceable[qid]})

        if valid_total <= 0:
            raise CommunityError("sample '" + sample["sample_id"] + "': no "
                                 "placeable reads; valid total must be positive")

        events: dict[int, list[tuple[float, float]]] = {}
        for eid, evs in raw_events.items():
            if not evs:
                continue
            merged: dict[float, float] = {}
            for pos, mass in evs:
                merged[pos] = merged.get(pos, 0.0) + mass / valid_total
            events[eid] = sorted(merged.items())
        measures.append(events)
        report.append({"sample_id": sample["sample_id"],
                       "job_id": sample["job_id"],
                       "excluded": excluded,
                       "valid_total": valid_total})

    return geo, measures, report


def _distal_mass(events_on_edge: list[tuple[float, float]], strict_after: float
                 ) -> float:
    """Mass at positions strictly farther from the parent than strict_after."""
    total = 0.0
    for pos, mass in events_on_edge:
        if pos > strict_after:
            total += mass
    return total


def pair_distance(geo: dict, m1: dict, m2: dict) -> dict:
    """KR (p=1) distance with per-edge non-negative contributions.

    On edge e the distal subtree mass at coordinate x (0 at parent, L at
    child) equals the full child-side subtree mass plus the mass sitting on e
    strictly distal of x.  The contribution is the integral over x of the
    absolute difference of the two samples' distal masses; events make that
    difference piecewise constant.
    """
    subtree = geo["subtree"]
    lengths = geo["lengths"]

    def edge_total(measure, eid):
        return sum(m for _, m in measure.get(eid, []))

    d1: dict[int, float] = {}
    d2: dict[int, float] = {}
    for eid in geo["edge_order"]:
        child_total_1 = sum(d1.get(c, 0.0) for c in subtree[eid])
        child_total_2 = sum(d2.get(c, 0.0) for c in subtree[eid])
        d1[eid] = child_total_1 + edge_total(m1, eid)
        d2[eid] = child_total_2 + edge_total(m2, eid)

    contributions: dict[int, float] = {}
    distance = 0.0
    for eid, length in lengths.items():
        ev1 = m1.get(eid, [])
        ev2 = m2.get(eid, [])
        breakpoints = sorted({0.0, length} | {p for p, _ in ev1}
                             | {p for p, _ in ev2})
        child_mass_1 = d1[eid] - edge_total(m1, eid)
        child_mass_2 = d2[eid] - edge_total(m2, eid)
        edge_contrib = 0.0
        for a, b in zip(breakpoints, breakpoints[1:]):
            if b <= a:
                continue
            distal_1 = child_mass_1 + _distal_mass(ev1, a)
            distal_2 = child_mass_2 + _distal_mass(ev2, a)
            edge_contrib += (b - a) * abs(distal_1 - distal_2)
        edge_contrib = max(edge_contrib, 0.0)
        contributions[eid] = edge_contrib
        distance += edge_contrib

    return {"distance": distance, "contributions": contributions}


def distance_matrix(geo: dict, measures: list[dict]) -> list[list[float]]:
    n = len(measures)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = pair_distance(geo, measures[i], measures[j])["distance"]
            matrix[i][j] = matrix[j][i] = d
    return matrix
