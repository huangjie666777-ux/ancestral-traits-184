"""Self-test: run a placement in-process and sanity-check the results."""

import json
import pathlib

from fastapi.testclient import TestClient

from app.tree import parse_annotated_newick, quote_name
from app.main import app

ROOT = pathlib.Path(__file__).resolve().parent.parent
EX = ROOT / "examples"


def main() -> None:
    client = TestClient(app)
    payload = {
        "reference_fasta": (EX / "reference.fasta").read_text(),
        "query_fasta": (EX / "queries.fasta").read_text(),
        "newick": (EX / "tree.nwk").read_text(),
    }
    resp = client.post("/place", json=payload)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert "q3" in data["unplaceable"], data
    for qid in ("q1", "q2"):
        places = data["placements"][qid]
        assert len(places) == data["edge_count"] == 7
        weights = [p["like_weight_ratio"] for p in places]
        assert abs(sum(weights) - 1.0) < 1e-9, weights
        lls = [p["likelihood"] for p in places]
        assert lls == sorted(lls, reverse=True)
        for p in places:
            assert 0.0 <= p["pendant_length"] <= 2.0
            assert p["distal_length"] >= 0.0
    # q1 matches refA exactly: its best edge should touch refA's pendant edge
    best_q1 = data["placements"]["q1"][0]
    print("q1 best edge:", best_q1)

    dl = client.get("/place/" + data["job_id"] + "/jplace")
    assert dl.status_code == 200
    doc = json.loads(dl.content)
    assert doc["version"] == 3
    assert doc["fields"] == ["edge_num", "likelihood", "like_weight_ratio",
                             "distal_length", "pendant_length"]
    assert "{" in doc["tree"] and "}" in doc["tree"]
    assert len(doc["placements"]) == 2

    # rejection path: illegal character must be located and reject the batch
    bad = dict(payload)
    bad["query_fasta"] = ">qX\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGX\n"
    resp = client.post("/place", json=bad)
    assert resp.status_code == 422
    problems = resp.json()["detail"]["problems"]
    assert any("qX" in p and "column 40" in p for p in problems), problems
    print("rejection OK:", problems[0])

    # RNA uracil must not be silently accepted as DNA
    bad_u = dict(payload)
    bad_u["query_fasta"] = ">qU\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGU\n"
    resp = client.post("/place", json=bad_u)
    assert resp.status_code == 422
    problems_u = resp.json()["detail"]["problems"]
    assert any("qU" in p and "column 40" in p and "'U'" in p
               for p in problems_u), problems_u
    print("U rejection OK:", problems_u[0])

    # ---- community KR distances on the stored job ----
    tree_text = doc["tree"]
    geo = parse_annotated_newick(tree_text)
    assert sorted(geo["lengths"]) == list(range(data["edge_count"]))
    assert all(v > 0 for v in geo["lengths"].values())

    # Newick labels containing commas/parens/quotes must survive quoting
    weird = ("(('we,ird':0.1{0},'x''y':0.12{1}):0.08{2},"
             "('r c':0.11{3},d:0.09{4}):0.07{5},e:0.2{6});")
    g2 = parse_annotated_newick(weird)
    assert g2["leaf"][0] == "we,ird", g2["leaf"]
    assert g2["leaf"][1] == "x'y", g2["leaf"]
    assert g2["lengths"][6] == 0.2
    assert quote_name("a'b,c") == "'a''b,c'"

    job_id = data["job_id"]
    comp = {
        "samples": [
            {"sample_id": "site1", "job_id": job_id,
             "counts": {"q1": 8, "q3": 3}},
            {"sample_id": "site2", "job_id": job_id,
             "counts": {"q1": 2, "q2": 6, "q3": 1}},
            {"sample_id": "site3", "job_id": job_id,
             "counts": {"q1": 5, "q2": 3}},
        ]
    }
    resp = client.post("/community/distances", json=comp)
    assert resp.status_code == 200, resp.text
    cdata = resp.json()
    assert cdata["sample_ids"] == ["site1", "site2", "site3"]
    mat = cdata["matrix"]
    assert len(mat) == 3 and all(len(row) == 3 for row in mat)
    for i in range(3):
        assert mat[i][i] == 0.0
        for j in range(3):
            assert mat[i][j] == mat[j][i]
            assert mat[i][j] >= 0.0
    assert mat[0][1] > 0.0 and mat[0][2] > 0.0
    for rep in cdata["samples"]:
        assert rep["valid_total"] > 0
    assert cdata["samples"][0]["excluded"] == [
        {"id": "q3", "count": 3,
         "reason": data["unplaceable"]["q3"]}]
    for pair in cdata["pairs"]:
        i, j = pair["i"], pair["j"]
        edges = [c["edge_num"] for c in pair["edge_contributions"]]
        assert edges == list(range(data["edge_count"]))
        contribs = [c["contribution"] for c in pair["edge_contributions"]]
        assert all(c >= 0.0 for c in contribs)
        assert abs(sum(contribs) - pair["distance"]) < 1e-9
        assert abs(pair["distance"] - mat[i][j]) < 1e-12
    print("community matrix:", [[round(v, 6) for v in row] for row in mat])

    # identical compositions (up to scaling) -> distance 0
    same = {"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 3, "q2": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 6, "q2": 2}}]}
    resp = client.post("/community/distances", json=same)
    assert resp.status_code == 200, resp.text
    assert resp.json()["matrix"][0][1] == 0.0

    def reject(body, fragment):
        r = client.post("/community/distances", json=body)
        assert r.status_code == 422, r.text
        assert fragment in r.json()["detail"]["problem"], r.text
        print("community rejected:", r.json()["detail"]["problem"])

    reject({"samples": [{"sample_id": "a", "job_id": job_id,
                         "counts": {"q1": 1}}]}, "2..10 samples")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}}]},
        "duplicate sample_id")
    reject({"samples": [
        {"sample_id": "a", "job_id": "nope", "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "unknown job_id")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"qX": 1}}]},
        "unknown query id")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": -1}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "non-negative integer")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1.5}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "non-negative integer")
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q3": 2}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "valid total must be positive")

    # a second job on a different reference tree must be rejected pairwise
    other = dict(payload)
    other["newick"] = ("((refA:0.1,refB:0.12):0.08,"
                       "(refC:0.11,refD:0.09):0.07,refE:0.31);")
    r2 = client.post("/place", json=other)
    assert r2.status_code == 200, r2.text
    job2 = r2.json()["job_id"]
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"q1": 1}},
        {"sample_id": "b", "job_id": job2, "counts": {"q1": 1}}]},
        "different annotated reference tree")
    # ---- deep-tree distal mass must be accumulated exactly once ----
    from app.community import pair_distance
    deep = ("((a:0.1{0},b:0.1{1}):0.2{2},"
            "(c:0.1{3},(d:0.1{4},e:0.1{5}):0.2{6}):0.3{7});")
    gdeep = parse_annotated_newick(deep)
    assert gdeep["subtree"][6] == [4, 5], gdeep["subtree"]
    assert gdeep["subtree"][7] == [3, 6], gdeep["subtree"]
    m_d = {4: [(0.1, 1.0)]}   # all mass at the tip of leaf d
    m_e = {5: [(0.1, 1.0)]}   # all mass at the tip of leaf e
    d_de = pair_distance(gdeep, m_d, m_e)["distance"]
    assert abs(d_de - 0.2) < 1e-12, d_de  # 0.1 up + 0.1 down, nothing more
    m_a = {0: [(0.1, 1.0)]}
    d_ad = pair_distance(gdeep, m_a, m_d)["distance"]
    # a up 0.1 + edge2 0.2 + edge7 0.3 + edge6 0.2 + d down 0.1 = 0.9
    assert abs(d_ad - 0.9) < 1e-12, d_ad
    print("deep-tree distal mass OK:", round(d_de, 6), round(d_ad, 6))

    # zero-count entries with unknown query ids must not slip through
    reject({"samples": [
        {"sample_id": "a", "job_id": job_id, "counts": {"qX": 0}},
        {"sample_id": "b", "job_id": job_id, "counts": {"q1": 1}}]},
        "unknown query id")

    # ---- batch-constrained PERMANOVA ----
    def pm_samples(groups, batches=None):
        comps = [{"q1": 8, "q2": 1}, {"q1": 7, "q2": 2},
                 {"q1": 1, "q2": 8}, {"q1": 2, "q2": 7}]
        out = []
        for i, (counts, g) in enumerate(zip(comps, groups)):
            s = {"sample_id": "s" + str(i), "job_id": job_id,
                 "counts": counts, "group": g}
            if batches is not None:
                s["batch"] = batches[i]
            out.append(s)
        return out

    # 4 samples, 2 groups of 2 -> 6 distinct assignments <= 999 -> exact
    body = {"samples": pm_samples(["T", "T", "C", "C"]),
            "permutations": 999, "seed": 7}
    resp = client.post("/community/permanova", json=body)
    assert resp.status_code == 200, resp.text
    pm = resp.json()
    assert pm["sample_ids"] == ["s0", "s1", "s2", "s3"]
    assert pm["groups"] == ["T", "T", "C", "C"]
    assert pm["batches"] is None
    assert pm["method"] == "exact" and pm["permutations"] == 6
    assert pm["n_assignments"] == 6
    assert pm["df_between"] == 1 and pm["df_within"] == 2
    assert abs(pm["ss_total"] - (pm["ss_within"] + pm["ss_between"])) < 1e-12
    assert 0.0 < pm["R2"] <= 1.0
    exp_f = ((pm["ss_between"] / pm["df_between"])
             / (pm["ss_within"] / pm["df_within"]))
    assert abs(pm["F"] - exp_f) < 1e-12
    assert abs(pm["R2"] - pm["ss_between"] / pm["ss_total"]) < 1e-12
    k = round(pm["p_value"] * 6)
    assert abs(pm["p_value"] - k / 6) < 1e-12 and 1 <= k <= 6
    print("permanova exact:", {k2: pm[k2] for k2 in
          ("F", "R2", "p_value", "method", "permutations")})

    # same request + same seed reproduces exactly (random path)
    body_r = {"samples": pm_samples(["T", "T", "C", "C"]),
              "permutations": 5, "seed": 42}
    r1 = client.post("/community/permanova", json=body_r).json()
    r2 = client.post("/community/permanova", json=body_r).json()
    assert r1 == r2 and r1["method"] == "random"
    assert r1["permutations"] == 5
    num = round(r1["p_value"] * 6)
    assert abs(r1["p_value"] - num / 6) < 1e-12 and 1 <= num <= 6
    print("permanova random reproducible:", r1["p_value"])

    # batches: labels only swap within batches; batch counts preserved
    body_b = {"samples": pm_samples(["T", "C", "T", "C"],
                                    ["b1", "b1", "b2", "b2"]),
              "permutations": 999, "seed": 1}
    rb = client.post("/community/permanova", json=body_b)
    assert rb.status_code == 200, rb.text
    pmb = rb.json()
    assert pmb["batches"] == ["b1", "b1", "b2", "b2"]
    assert pmb["n_assignments"] == 4 and pmb["permutations"] == 4
    print("permanova batched exact:", pmb["p_value"])

    def pm_reject(body, fragment):
        r = client.post("/community/permanova", json=body)
        assert r.status_code == 422, r.text
        assert fragment in r.json()["detail"]["problem"], r.text
        print("permanova rejected:", r.json()["detail"]["problem"])

    pm_reject({"samples": pm_samples(["T", "T", "C"])[:3],
               "permutations": 99, "seed": 1}, "4..10 samples")
    pm_reject({"samples": pm_samples(["T", "T", "T", "C"]),
               "permutations": 99, "seed": 1}, "at least 2")
    pm_reject({"samples": pm_samples(["T", "T", "T", "T"]),
               "permutations": 99, "seed": 1}, "at least 2 groups")
    pm_reject({"samples": pm_samples(["T", "T", "C", "C"]),
               "permutations": 0, "seed": 1}, "permutations")
    pm_reject({"samples": pm_samples(["T", "T", "C", "C"]),
               "permutations": 10000, "seed": 1}, "permutations")
    pm_reject({"samples": pm_samples(["T", "T", "C", "C"]),
               "permutations": 99, "seed": "x"}, "seed")
    pm_reject({"samples": pm_samples(["T", "T", "C", "C"],
                                     ["b1", "b1", "b2", None]),
               "permutations": 99, "seed": 1}, "every sample must")
    # each batch internally uniform -> no legal exchange exists
    pm_reject({"samples": pm_samples(["T", "T", "C", "C"],
                                     ["b1", "b1", "b2", "b2"]),
               "permutations": 99, "seed": 1}, "no legal label exchange")
    # identical compositions everywhere -> total SS is zero
    pm_reject({"samples": [
        {"sample_id": "u" + str(i), "job_id": job_id,
         "counts": {"q1": 3, "q2": 1}, "group": g}
        for i, g in enumerate(["T", "T", "C", "C"])],
        "permutations": 99, "seed": 1}, "total sum of squares is zero")
    # groups internally identical but distinct -> within SS is zero
    pm_reject({"samples": [
        {"sample_id": "w0", "job_id": job_id,
         "counts": {"q1": 3, "q2": 1}, "group": "T"},
        {"sample_id": "w1", "job_id": job_id,
         "counts": {"q1": 6, "q2": 2}, "group": "T"},
        {"sample_id": "w2", "job_id": job_id,
         "counts": {"q1": 1, "q2": 4}, "group": "C"},
        {"sample_id": "w3", "job_id": job_id,
         "counts": {"q1": 2, "q2": 8}, "group": "C"}],
        "permutations": 99, "seed": 1}, "within-group sum of squares is zero")

    # ---- ancestral state reconstruction on the stored job ----
    from app.ancestral import _build_nodes, _reroot
    import itertools

    anc = {
        "job_id": job_id,
        "states": ["tolerant", "sensitive"],
        "leaf_states": {
            "refA": ["tolerant"],
            "refB": ["tolerant"],
            "refrefC": ["sensitive"],
            "refD": ["tolerant", "sensitive"],
            "refE": ["sensitive"],
        },
        "outgroup": "refE",
        "transition_costs": [[0, 2], [3, 0]],
    }
    anc["leaf_states"] = {
        "refA": ["tolerant"], "refB": ["tolerant"], "refC": ["sensitive"],
        "refD": ["tolerant", "sensitive"], "refE": ["sensitive"]}
    resp = client.post("/ancestral/reconstruct", json=anc)
    assert resp.status_code == 200, resp.text
    ar = resp.json()
    assert ar["feasible"] is True
    assert isinstance(ar["optimal_histories"], str)
    assert ar["min_cost"] >= 0

    # brute-force cross-check over every complete assignment
    geo_a = parse_annotated_newick(tree_text)
    nodes_a = _build_nodes(geo_a)
    root_a, children_a, parent_edge_a, endpoints_a = _reroot(nodes_a, "refE")
    n_nodes = len(nodes_a)
    allowed_idx = []
    for nd in nodes_a:
        if nd["leaf"] is None:
            allowed_idx.append([0, 1])
        else:
            allowed_idx.append([ar["states"].index(s)
                                for s in anc["leaf_states"][nd["leaf"]]])
    mat = anc["transition_costs"]
    best_cost = None
    best_assigns = []
    for combo in itertools.product(*allowed_idx):
        total = 0
        ok = True
        for e, (u, v) in endpoints_a.items():
            w = mat[combo[u]][combo[v]]
            if w is None:
                ok = False
                break
            total += w
        if not ok:
            continue
        if best_cost is None or total < best_cost:
            best_cost, best_assigns = total, [combo]
        elif total == best_cost:
            best_assigns.append(combo)
    assert ar["min_cost"] == best_cost, (ar["min_cost"], best_cost)
    assert int(ar["optimal_histories"]) == len(best_assigns), (
        ar["optimal_histories"], len(best_assigns))
    node_by_id = {}
    for i, nd in enumerate(nodes_a):
        nid = ("leaf:" + nd["leaf"] if nd["leaf"] is not None
               else "internal:" + ",".join(str(e) for e in nd["edges"]))
        node_by_id[nid] = i
    for rep in ar["nodes"]:
        i = node_by_id[rep["node_id"]]
        expect = sorted({ar["states"][c[i]] for c in best_assigns},
                        key=ar["states"].index)
        assert rep["possible_states"] == expect, (rep, expect)
    for rep in ar["edges"]:
        u, v = endpoints_a[rep["edge_num"]]
        expect = sorted({(ar["states"][c[u]], ar["states"][c[v]])
                         for c in best_assigns})
        got = sorted(tuple(p) for p in rep["possible_pairs"])
        assert got == expect, (rep, expect)
        kinds = {a != b for a, b in expect}
        want = ("always_change" if kinds == {True} else
                "never_change" if kinds == {False} else "may_change")
        assert rep["change"] == want, (rep, want)
    root_entries = [n for n in ar["nodes"] if n["is_root"]]
    assert len(root_entries) == 1
    assert root_entries[0]["parent"] is None
    assert root_entries[0]["parent_edge"] is None
    assert ar["root"] == root_entries[0]["node_id"]
    assert all(n["parent"] is not None for n in ar["nodes"]
               if not n["is_root"])
    leaf_entries = [n for n in ar["nodes"] if n["kind"] == "leaf"]
    assert sorted(n["leaf_id"] for n in leaf_entries) == \
        ["refA", "refB", "refC", "refD", "refE"]
    assert sorted(e["edge_num"] for e in ar["edges"]) == list(range(7))
    print("ancestral:", ar["min_cost"], "cost,",
          ar["optimal_histories"], "optimal histories, root", ar["root"])

    # unknown trait (all states) inflates the count of optimal histories
    anc_unknown = dict(anc, leaf_states={
        "refA": ["tolerant"], "refB": ["tolerant"], "refC": ["sensitive"],
        "refD": ["tolerant", "sensitive"],
        "refE": ["tolerant", "sensitive"]})
    resp = client.post("/ancestral/reconstruct", json=anc_unknown)
    assert resp.status_code == 200, resp.text
    au = resp.json()
    assert au["min_cost"] <= ar["min_cost"] + 3
    assert int(au["optimal_histories"]) >= 1
    print("ancestral with unknown outgroup trait:", au["optimal_histories"])

    # forbidden transition can make every history infeasible
    anc_null = dict(anc, transition_costs=[[0, None], [None, 0]])
    resp = client.post("/ancestral/reconstruct", json=anc_null)
    assert resp.status_code == 200, resp.text
    inf = resp.json()
    assert inf["feasible"] is False
    assert "no feasible history" in inf["reason"]
    assert "nodes" not in inf and "edges" not in inf
    print("infeasible reason:", inf["reason"])

    def anc_reject(body, fragment):
        r = client.post("/ancestral/reconstruct", json=body)
        assert r.status_code == 422, r.text
        assert fragment in r.json()["detail"]["problem"], r.text
        print("ancestral rejected:", r.json()["detail"]["problem"])

    anc_reject(dict(anc, job_id="nope"), "unknown job_id")
    anc_reject(dict(anc, states=["tolerant"]), "2..5 states")
    anc_reject(dict(anc, states=["a", "a"]), "duplicate state name")
    anc_reject(dict(anc, leaf_states={k: v for k, v in
                                      anc["leaf_states"].items()
                                      if k != "refB"}),
               "leaf 'refB' has no allowed-state set")
    bad_leaf = dict(anc["leaf_states"], refZ=["tolerant"])
    anc_reject(dict(anc, leaf_states=bad_leaf), "unknown leaf 'refZ'")
    bad_state = dict(anc["leaf_states"], refA=["tolerant", "weird"])
    anc_reject(dict(anc, leaf_states=bad_state), "unknown state 'weird'")
    dup_state = dict(anc["leaf_states"],
                     refA=["tolerant", "tolerant"])
    anc_reject(dict(anc, leaf_states=dup_state), "duplicate state 'tolerant'")
    empty_set = dict(anc["leaf_states"], refA=[])
    anc_reject(dict(anc, leaf_states=empty_set), "non-empty list")
    anc_reject(dict(anc, outgroup="refZ"), "unknown outgroup leaf 'refZ'")
    anc_reject(dict(anc, transition_costs=[[0, 1], [1]]), "list of 2 entries")
    anc_reject(dict(anc, transition_costs=[[0, 1], [1, 1]]),
               "the diagonal must be 0")
    anc_reject(dict(anc, transition_costs=[[0, True], [1, 0]]),
               "booleans are not allowed")
    anc_reject(dict(anc, transition_costs=[[0, 1.5], [1, 0]]),
               "non-negative integer or null")
    anc_reject(dict(anc, transition_costs=[[0, -1], [1, 0]]),
               "non-negative integer or null")
    anc_reject(dict(anc, transition_costs=[[None, 1], [1, 0]]),
               "the diagonal must be 0, not null")

    # stored job must remain untouched by reconstruction
    dl2 = client.get("/place/" + job_id + "/jplace")
    assert dl2.status_code == 200
    assert json.loads(dl2.content) == doc

    print("SELFTEST PASSED")


if __name__ == "__main__":
    main()
