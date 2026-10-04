"""Batch-constrained PERMANOVA on the community KR distance matrix.

The design (per-sample group, optional batch) and the permutation scheme
live here; placement jobs, mass construction and the KR distance matrix
are reused from app.community unchanged (no re-placement, no mutation).

Sums of squares follow the distance-based (Anderson) definition:
  SS_total  = sum_{i<j} d_ij^2 / n
  SS_within = sum_g (sum_{i<j in g} d_ij^2 / n_g)
  SS_between = SS_total - SS_within
pseudo-F = (SS_between / (g-1)) / (SS_within / (n-g));
R2 = SS_between / SS_total.

Only group labels are permuted; distances and group sizes stay fixed.
With batches, labels are shuffled within each batch only, preserving
per-batch group counts.  If the number of distinct legal assignments
does not exceed the requested permutations, all of them (observed
included) are enumerated and p is the fraction with F >= F_obs;
otherwise the requested number of seeded random draws (repeats allowed) give
p = (extreme + 1) / (permutations + 1).
"""

from __future__ import annotations

import math
import random
from itertools import product

from .community import (CommunityError, build_measures, distance_matrix,
                        resolve_jobs, validate_request)

MIN_SAMPLES = 4
MAX_SAMPLES = 10
MIN_PERMUTATIONS = 1
MAX_PERMUTATIONS = 9999
_EPS = 1e-9


def _validate_design(payload) -> tuple:
    """Validate samples + group/batch design + permutation settings."""
    if not isinstance(payload, dict) or "samples" not in payload:
        raise CommunityError("request body must contain 'samples'")
    raw = payload["samples"]
    if not isinstance(raw, list):
        raise CommunityError("'samples' must be a list")
    if not MIN_SAMPLES <= len(raw) <= MAX_SAMPLES:
        raise CommunityError("expected " + str(MIN_SAMPLES) + ".." +
                             str(MAX_SAMPLES) + " samples, got " + str(len(raw)))
    samples = validate_request(payload)  # sample_id/job_id/counts checks

    any_batch = any("batch" in s for s in samples)
    groups: list = []
    batches: list = []
    for idx, sample in enumerate(samples):
        where = "sample[" + str(idx) + "] ('" + sample["sample_id"] + "')"
        group = sample.get("group")
        if not isinstance(group, str) or not group:
            raise CommunityError(where + ".group must be a non-empty string")
        groups.append(group)
        if any_batch:
            batch = sample.get("batch")
            if not isinstance(batch, str) or not batch:
                raise CommunityError(
                    where + ".batch must be a non-empty string: when any "
                    "sample gives a batch, every sample must")
            batches.append(batch)
        else:
            batches.append(None)

    sizes: dict = {}
    for g in groups:
        sizes[g] = sizes.get(g, 0) + 1
    if len(sizes) < 2:
        raise CommunityError("need at least 2 groups, got " + str(len(sizes)))
    for g in sorted(sizes):
        if sizes[g] < 2:
            raise CommunityError("group '" + g + "' has only " +
                                 str(sizes[g]) + " sample(s); every group "
                                 "needs at least 2")

    perm = payload.get("permutations")
    if (not isinstance(perm, int) or isinstance(perm, bool)
            or not MIN_PERMUTATIONS <= perm <= MAX_PERMUTATIONS):
        raise CommunityError("'permutations' must be an integer in " +
                             str(MIN_PERMUTATIONS) + ".." +
                             str(MAX_PERMUTATIONS) + ", got " + repr(perm))
    seed = payload.get("seed")
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise CommunityError("'seed' must be an integer, got " + repr(seed))
    return samples, groups, batches, perm, seed


def _ss_total(d2: list) -> float:
    n = len(d2)
    return sum(d2[i][j] for i in range(n) for j in range(i + 1, n)) / n


def _ss_within(labels: list, d2: list) -> float:
    members: dict = {}
    for i, g in enumerate(labels):
        members.setdefault(g, []).append(i)
    total = 0.0
    for idxs in members.values():
        pair_sum = 0.0
        for a in range(len(idxs)):
            for b in range(a + 1, len(idxs)):
                pair_sum += d2[idxs[a]][idxs[b]]
        total += pair_sum / len(idxs)
    return total


def _pseudo_f(labels: list, d2: list, ss_total: float,
              df_between: int, df_within: int) -> float:
    ss_w = _ss_within(labels, d2)
    if ss_w <= 0.0:
        return math.inf
    ss_b = ss_total - ss_w
    return (ss_b / df_between) / (ss_w / df_within)


def _unique_permutations(seq: list):
    """Yield each distinct permutation of a multiset exactly once."""
    items = sorted(seq)

    def rec(prefix, rest):
        if not rest:
            yield prefix
            return
        for i in range(len(rest)):
            if i > 0 and rest[i] == rest[i - 1]:
                continue
            yield from rec(prefix + [rest[i]], rest[:i] + rest[i + 1:])

    yield from rec([], items)


def _batch_blocks(batches: list) -> list:
    blocks: dict = {}
    for i, b in enumerate(batches):
        blocks.setdefault(b, []).append(i)
    return list(blocks.values())


def _n_assignments(groups: list, blocks: list) -> int:
    """Distinct legal assignments: product of per-batch multinomials."""
    total = 1
    for idxs in blocks:
        counts: dict = {}
        for i in idxs:
            counts[groups[i]] = counts.get(groups[i], 0) + 1
        ways = math.factorial(len(idxs))
        for c in counts.values():
            ways //= math.factorial(c)
        total *= ways
    return total


def run_permanova(payload: dict, jobs: dict) -> dict:
    samples, groups, batches, perm, seed = _validate_design(payload)
    resolved = resolve_jobs(samples, jobs)
    geo, measures, _report = build_measures(samples, resolved)
    matrix = distance_matrix(geo, measures)
    d2 = [[v * v for v in row] for row in matrix]
    n = len(samples)

    ss_total = _ss_total(d2)
    if ss_total <= 0.0:
        raise CommunityError(
            "total sum of squares is zero: all pairwise KR distances "
            "vanish, so no community difference can be tested")
    ss_within = _ss_within(groups, d2)
    if ss_within <= 0.0:
        raise CommunityError(
            "within-group sum of squares is zero: pseudo-F is undefined "
            "(every group is internally identical)")

    n_groups = len(set(groups))
    df_between = n_groups - 1
    df_within = n - n_groups
    ss_between = ss_total - ss_within
    f_obs = (ss_between / df_between) / (ss_within / df_within)
    r2 = ss_between / ss_total

    blocks = _batch_blocks(batches)
    n_assign = _n_assignments(groups, blocks)
    if n_assign < 2:
        raise CommunityError(
            "no legal label exchange changes the grouping (within every "
            "batch each group label is uniform); a permutation test is "
            "impossible and cross-batch shuffling is not a valid test")

    if n_assign <= perm:
        method = "exact"
        extreme = 0
        used = 0
        per_block = [[groups[i] for i in idxs] for idxs in blocks]
        for combo in product(*[_unique_permutations(labels)
                               for labels in per_block]):
            labels = list(groups)
            for idxs, permuted in zip(blocks, combo):
                for pos, i in enumerate(idxs):
                    labels[i] = permuted[pos]
            f_val = _pseudo_f(labels, d2, ss_total, df_between, df_within)
            used += 1
            if f_val >= f_obs - _EPS:
                extreme += 1
        p_value = extreme / used
    else:
        method = "random"
        rng = random.Random(seed)
        extreme = 0
        used = perm
        for _ in range(perm):
            labels = list(groups)
            for idxs in blocks:
                permuted = [groups[i] for i in idxs]
                rng.shuffle(permuted)
                for pos, i in enumerate(idxs):
                    labels[i] = permuted[pos]
            f_val = _pseudo_f(labels, d2, ss_total, df_between, df_within)
            if f_val >= f_obs - _EPS:
                extreme += 1
        p_value = (extreme + 1) / (perm + 1)

    return {
        "sample_ids": [s["sample_id"] for s in samples],
        "groups": groups,
        "batches": batches if any(b is not None for b in batches) else None,
        "ss_total": ss_total,
        "ss_within": ss_within,
        "ss_between": ss_between,
        "df_between": df_between,
        "df_within": df_within,
        "F": f_obs,
        "R2": r2,
        "p_value": p_value,
        "method": method,
        "permutations": used,
        "n_assignments": n_assign,
    }
