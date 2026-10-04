"""Self-test: run a placement in-process and sanity-check the results."""

import json
import pathlib

from fastapi.testclient import TestClient

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
    print("SELFTEST PASSED")


if __name__ == "__main__":
    main()
