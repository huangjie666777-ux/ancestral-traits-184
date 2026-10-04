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
    print("SELFTEST PASSED")


if __name__ == "__main__":
    main()
