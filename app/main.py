"""FastAPI entrypoint: validate, place, and deliver jplace results."""

from __future__ import annotations

import json
import uuid

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.responses import Response
from pydantic import BaseModel

from .jplacefmt import build_jplace
from .likelihood import DirectedMessages, obs_vector
from .placement import has_base_evidence, place_query
from .community import (CommunityError, build_measures, distance_matrix,
                        pair_distance, resolve_jobs, validate_request)
from .tree import from_phylotree, orient_and_serialize, renumber_edges
from .validation import (MAX_QUERIES, MAX_REFS, MIN_QUERIES, MIN_REFS,
                         SubmissionError, cross_validate, parse_tree,
                         validate_alignment)

app = FastAPI(title="eDNA placement backend", version="1.0.0")

JOBS: dict[str, dict] = {}


class PlacementRequest(BaseModel):
    reference_fasta: str
    query_fasta: str
    newick: str


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
