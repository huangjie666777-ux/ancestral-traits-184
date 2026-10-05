"""FastAPI entrypoint: validate, place, and deliver jplace results."""

from __future__ import annotations

import json
import uuid

import numpy as np
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import Response
from pydantic import BaseModel

from .ancstates import AncestralError, reconstruct
from .jplacefmt import build_jplace
from .likelihood import DirectedMessages, obs_vector
from .placement import has_base_evidence, place_query
from .community import (CommunityError, build_measures, distance_matrix,
                        pair_distance, resolve_jobs, validate_request)
from .permanova import run_permanova
from .tree import (from_phylotree, orient_and_serialize,
                   parse_annotated_newick, renumber_edges)
from .validation import (MAX_QUERIES, MAX_REFS, MIN_QUERIES, MIN_REFS,
                         SubmissionError, cross_validate, parse_tree,
                         validate_alignment)

app = FastAPI(title="eDNA placement backend", version="1.0.0")

JOBS: dict[str, dict] = {}


class PlacementRequest(BaseModel):
    reference_fasta: str
    query_fasta: str
    newick: str


def _check_state_list(value, where: str) -> list[str]:
    """2..5 unique non-empty state names, in declared order."""
    if not isinstance(value, list):
        raise AncestralError(where + " must be a list of state names")
    if not 2 <= len(value) <= 5:
        raise AncestralError(where + ": expected 2..5 states, got " +
                             str(len(value)))
    seen: set[str] = set()
    for i, name in enumerate(value):
        if not isinstance(name, str) or not name:
            raise AncestralError(where + "[" + str(i) +
                                 "]: state name must be a non-empty string")
        if name in seen:
            raise AncestralError(where + ": duplicate state name '" +
                                 name + "'")
        seen.add(name)
    return list(value)


def _pairs_hook(pairs):
    """Keep each object's raw (key, value) pairs visible for duplicate checks.

    Leaf-state maps arrive as tuples tagged '__raw_pairs__' so a repeated
    leaf key (which a plain dict would silently collapse) can be located and
    reject the whole request.
    """
    result: dict = {}
    for k, v in pairs:
        result[k] = v
    if pairs and all(isinstance(k, str) for k, _ in pairs):
        result["__raw_pairs__"] = [(k, v) for k, v in pairs]
    return result
    return result


def validate_ancestral_body(payload: dict):
    """Locate every structural problem; the caller rejects the whole batch."""
    job_id = payload.get("job_id")
    if not isinstance(job_id, str) or not job_id:
        raise AncestralError("job_id must be a non-empty string")
    if job_id not in JOBS:
        raise AncestralError("unknown job_id '" + job_id + "'")

    states = _check_state_list(payload.get("states"), "states")

    outgroup = payload.get("outgroup_leaf_id")
    if not isinstance(outgroup, str) or not outgroup:
        raise AncestralError("outgroup_leaf_id must be a non-empty string")

    leaf_states = payload.get("leaf_states")
    if not isinstance(leaf_states, dict):
        raise AncestralError("leaf_states must be an object mapping each "
                             "reference leaf id to a list of allowed states")
    ls_pairs = leaf_states.get("__raw_pairs__", [])
    dup_leaves = {k for k, _ in ls_pairs
                  if sum(1 for kk, _ in ls_pairs if kk == k) > 1}
    if dup_leaves:
        raise AncestralError("leaf_states: duplicate leaf entry '" +
                             sorted(dup_leaves)[0] + "'")

    tree_text = JOBS[job_id]["jplace"]["tree"]
    try:
        geo = parse_annotated_newick(tree_text)
    except ValueError as exc:
        raise AncestralError("stored reference tree is unreadable: " +
                             str(exc))
    ref_leaves = {leaf for leaf in geo["leaf"].values() if leaf is not None}

    if outgroup not in ref_leaves:
        raise AncestralError("outgroup_leaf_id '" + outgroup +
                             "' is not a leaf of the reference tree")

    state_set = set(states)
    problems: list[str] = []
    allowed: dict[str, frozenset] = {}
    for leaf, vals in list(leaf_states.items()):
        if leaf == "__raw_pairs__":
            continue
        if not isinstance(leaf, str) or not leaf:
            problems.append("leaf_states: leaf id must be a non-empty string")
            continue
        if leaf not in ref_leaves:
            problems.append("leaf_states: leaf '" + leaf +
                            "' is not a leaf of the reference tree")
            continue
        if not isinstance(vals, list) or not vals:
            problems.append("leaf_states['" + leaf +
                            "']: allowed states must be a non-empty list")
            continue
        picked: list[str] = []
        seen: set[str] = set()
        for j, st in enumerate(vals):
            if not isinstance(st, str) or not st:
                problems.append("leaf_states['" + leaf + "'][" + str(j) +
                                "]: state must be a non-empty string")
            elif st not in state_set:
                problems.append("leaf_states['" + leaf + "']: unknown state "
                                "'" + st + "'")
            elif st in seen:
                problems.append("leaf_states['" + leaf +
                                "']: duplicate state '" + st + "'")
            else:
                seen.add(st)
                picked.append(st)
        if picked:
            allowed[leaf] = frozenset(states.index(s) for s in picked)
    missing = sorted(ref_leaves -
                     {k for k in leaf_states if k != "__raw_pairs__"})
    for leaf in missing:
        problems.append("leaf_states: missing entry for reference leaf '" +
                        leaf + "'")

    matrix = payload.get("cost_matrix")
    k = len(states)
    cost = None
    if not isinstance(matrix, list):
        problems.append("cost_matrix must be a list of rows")
    else:
        if len(matrix) != k:
            problems.append("cost_matrix: expected " + str(k) + " rows "
                            "(states order), got " + str(len(matrix)))
        cost = []
        for i, row in enumerate(matrix):
            if not isinstance(row, list):
                problems.append("cost_matrix[" + str(i) + "] must be a list")
                continue
            if len(row) != k:
                problems.append("cost_matrix[" + str(i) + "]: expected " +
                                str(k) + " columns, got " + str(len(row)))
                continue
            out_row = []
            for j, cell in enumerate(row):
                if cell is None:
                    out_row.append(None)
                elif isinstance(cell, bool) or not isinstance(cell, int) \
                        or cell < 0:
                    problems.append("cost_matrix[" + str(i) + "][" + str(j) +
                                    "]: cost must be a non-negative integer "
                                    "or null, got " + repr(cell))
                    out_row.append(None)
                else:
                    out_row.append(cell)
            if len(out_row) == k:
                if out_row[i] != 0:
                    problems.append("cost_matrix[" + str(i) + "][" + str(i) +
                                    "]: diagonal cost must be 0, got " +
                                    repr(out_row[i]))
            cost.append(out_row)

    if problems:
        raise AncestralError("; ".join(problems))
    return tree_text, states, allowed, cost


@app.post("/ancestral/states", status_code=200)
async def ancestral_states(request: Request):
    """Ancestral state reconstruction on a stored placement job's tree.

    Reads the job (never modified, placement never recomputed), roots the
    reference tree at the internal node adjacent to the outgroup leaf and
    minimises total Sankoff cost under leaf constraints.  Returns the minimum
    cost, the exact decimal-string count of all globally optimal complete
    assignments, the states possible at each node in *some* global optimum
    and the ordered parent->child state pairs possible on each edge.
    """
    raw = await request.body()
    try:
        payload = json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs_hook)
    except AncestralError:
        raise
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=422,
                            detail={"rejected": True,
                                    "problems": ["invalid JSON body: " +
                                                 str(exc)]})
    if not isinstance(payload, dict):
        raise HTTPException(status_code=422,
                            detail={"rejected": True,
                                    "problems": ["request body must be an "
                                                 "object"]})
    top_raw = payload.get("__raw_pairs__", [])
    dup_fields = {k for k, _ in top_raw
                  if not k.startswith("__") and
                  sum(1 for kk, _ in top_raw if kk == k) > 1}
    if dup_fields:
        raise HTTPException(status_code=422,
                            detail={"rejected": True,
                                    "problems": ["request body: duplicate "
                                                 "field '" +
                                                 sorted(dup_fields)[0] + "'"]})
    try:
        tree_text, states, allowed, cost = validate_ancestral_body(payload)
        result = reconstruct(tree_text, payload["outgroup_leaf_id"], states,
                             allowed, cost)
    except AncestralError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True,
                                    "problems": [str(exc)]})
    result["job_id"] = payload["job_id"]
    result["states"] = states
    result["outgroup_leaf_id"] = payload["outgroup_leaf_id"]
    return result


def run_placement(req: PlacementRequest) -> dict:
    problems: list[str] = []
    ref = validate_alignment(req.reference_fasta, "reference",
                             MIN_REFS, MAX_REFS, problems)
    qry = validate_alignment(req.query_fasta, "query",
                             MIN_QUERIES, MAX_QUERIES, problems)
    tree = parse_tree(req.newick, problems)
    if ref is not None and qry is not None and tree is not None:
        cross_validate(ref, qry, tree, problems)
    if problems:
        raise SubmissionError(problems)

    g = from_phylotree(tree)
    renumber_edges(g)
    tree_str = orient_and_serialize(g)

    name_to_node = {name: n for n, name in g.leaf_name.items()}
    leaf_obs = {}
    for sid, states in zip(ref.ids, ref.states):
        leaf_obs[name_to_node[sid]] = np.array([obs_vector(s) for s in states])
    dm = DirectedMessages(g, leaf_obs)

    placements_by_query: dict[str, list] = {}
    unplaceable: dict[str, str] = {}
    for sid, states in zip(qry.ids, qry.states):
        if not has_base_evidence(states):
            unplaceable[sid] = ("no valid base evidence: every column is a gap "
                                "or fully ambiguous, so the sequence carries no "
                                "phylogenetic signal and cannot be placed")
            continue
        query_obs = np.array([obs_vector(s) for s in states])
        placements_by_query[sid] = place_query(dm, g, query_obs)

    jplace = build_jplace(tree_str, placements_by_query)
    return {"jplace": jplace, "unplaceable": unplaceable,
            "edge_count": len(g.edges)}


@app.post("/place", status_code=200)
def place(req: PlacementRequest):
    try:
        result = run_placement(req)
    except SubmissionError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True, "problems": exc.problems})
    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = result
    summary = {
        "job_id": job_id,
        "jplace_url": "/place/" + job_id + "/jplace",
        "edge_count": result["edge_count"],
        "unplaceable": result["unplaceable"],
        "placements": {
            p["n"][0]: [
                dict(zip(result["jplace"]["fields"], row)) for row in p["p"]
            ]
            for p in result["jplace"]["placements"]
        },
    }
    return summary


@app.get("/place/{job_id}/jplace")
def download_jplace(job_id: str):
    if job_id not in JOBS:
        raise HTTPException(status_code=404, detail="unknown job_id")
    payload = json.dumps(JOBS[job_id]["jplace"], indent=2)
    return Response(
        content=payload,
        media_type="application/json",
        headers={"Content-Disposition":
                 "attachment; filename=placements_" + job_id + ".jplace"},
    )


@app.post("/community/distances", status_code=200)
def community_distances(payload: dict):
    """KR (p=1) distances between 2..10 read-count samples on one reference.

    Stored placement jobs are never modified or recomputed: placement mass is
    derived from each job's existing jplace (all candidate edges weighted by
    like_weight_ratio).  Any located problem (unknown job/query, illegal
    count, mismatched reference tree, non-positive valid total) rejects the
    whole batch with HTTP 422.
    """
    try:
        samples = validate_request(payload)
        jobs = resolve_jobs(samples, JOBS)
        geo, measures, report = build_measures(samples, jobs)
        matrix = distance_matrix(geo, measures)
        n = len(samples)
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                detail = pair_distance(geo, measures[i], measures[j])
                contrib_sum = sum(detail["contributions"].values())
                pairs.append({
                    "i": i, "j": j,
                    "sample_a": samples[i]["sample_id"],
                    "sample_b": samples[j]["sample_id"],
                    "distance": detail["distance"],
                    "contributions_sum": contrib_sum,
                    "edge_contributions": [
                        {"edge_num": eid,
                         "contribution": detail["contributions"][eid]}
                        for eid in sorted(geo["lengths"])
                    ],
                })
    except CommunityError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True, "problem": str(exc)})
    return {
        "sample_ids": [s["sample_id"] for s in samples],
        "tree": jobs[0]["jplace"]["tree"],
        "matrix": matrix,
        "pairs": pairs,
        "samples": report,
    }


@app.post("/community/permanova", status_code=200)
def community_permanova(payload: dict):
    """Batch-constrained PERMANOVA on KR distances between 4..10 samples.

    Reuses stored jobs and the community KR matrix; only group labels are
    permuted (within batches when given).  Illegal designs, degenerate
    sums of squares and designs with no legal label exchange are located
    and rejected with HTTP 422.
    """
    try:
        return run_permanova(payload, JOBS)
    except CommunityError as exc:
        raise HTTPException(status_code=422,
                            detail={"rejected": True, "problem": str(exc)})
