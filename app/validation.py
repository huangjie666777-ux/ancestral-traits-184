"""Input parsing and validation for reference/query FASTA and Newick tree.

All problems are collected with their location (record id, column, or tree
element) and the whole submission is rejected if any problem is found.
"""

from __future__ import annotations

from dataclasses import dataclass
from io import StringIO

from Bio import Phylo

# IUPAC ambiguity codes -> compatible nucleotide states (A, C, G, T).
IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "U": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"),
    "W": frozenset("AT"), "K": frozenset("GT"), "M": frozenset("AC"),
    "B": frozenset("CGT"), "D": frozenset("AGT"), "H": frozenset("ACT"),
    "V": frozenset("ACG"), "N": frozenset("ACGT"),
    "-": frozenset(), ".": frozenset(), "?": frozenset(),
}

MIN_REFS, MAX_REFS = 3, 20
MIN_QUERIES, MAX_QUERIES = 1, 5
MAX_COLUMNS = 500


class SubmissionError(Exception):
    """Raised with located problems; the whole submission is rejected."""

    def __init__(self, problems: list[str]):
        self.problems = problems
        super().__init__("; ".join(problems))


@dataclass
class Alignment:
    ids: list[str]
    states: list[list[frozenset]]  # states[i][col]; empty set = gap, no evidence
    n_cols: int


def parse_fasta(text: str, label: str, problems: list[str]):
    ids: list[str] = []
    seqs: list[str] = []
    current_id: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal current_id, chunks
        if current_id is not None:
            ids.append(current_id)
            seqs.append("".join(chunks))
        current_id, chunks = None, []

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if not header:
                problems.append(label + " FASTA line " + str(lineno) + ": empty sequence id")
                current_id = "<line " + str(lineno) + ">"
            else:
                current_id = header.split()[0]
        elif current_id is None:
            problems.append(label + " FASTA line " + str(lineno) +
                            ": sequence data before any header")
        else:
            chunks.append(line)
    flush()
    return ids, seqs


def validate_alignment(text, label, min_n, max_n, problems):
    ids, seqs = parse_fasta(text, label, problems)
    if not ids:
        problems.append(label + " FASTA: no sequences found")
        return None
    if not (min_n <= len(ids) <= max_n):
        problems.append(label + " FASTA: expected " + str(min_n) + ".." + str(max_n) +
                        " sequences, got " + str(len(ids)))
    seen = set()
    for sid in ids:
        if sid in seen:
            problems.append(label + " FASTA: duplicate sequence id '" + sid + "'")
        seen.add(sid)
    n_cols = len(seqs[0]) if seqs else 0
    if n_cols == 0:
        problems.append(label + " FASTA record '" + ids[0] + "': empty sequence")
    if n_cols > MAX_COLUMNS:
        problems.append(label + " FASTA: alignment has " + str(n_cols) +
                        " columns, maximum is " + str(MAX_COLUMNS))
    states = []
    for sid, seq in zip(ids, seqs):
        if len(seq) == 0:
            problems.append(label + " FASTA record '" + sid + "': empty sequence")
            continue
        if len(seq) != n_cols:
            problems.append(label + " FASTA record '" + sid + "': length " + str(len(seq)) +
                            " != " + str(n_cols) + " (all sequences must be aligned "
                            "to equal length)")
        row = []
        for col, ch in enumerate(seq.upper(), start=1):
            if ch not in IUPAC:
                problems.append(label + " FASTA record '" + sid + "' column " + str(col) +
                                ": illegal character '" + ch + "'")
                row.append(frozenset())
            else:
                row.append(IUPAC[ch])
        states.append(row)
    if problems:
        return None
    return Alignment(ids=ids, states=states, n_cols=n_cols)


def parse_tree(newick: str, problems: list[str]):
    """Parse the Newick tree; validate binary shape and branch lengths."""
    try:
        tree = Phylo.read(StringIO(newick.strip()), "newick")
    except Exception as exc:  # noqa: BLE001 - surface any parse failure
        problems.append("Newick tree: parse error: " + str(exc))
        return None
    clades = list(tree.find_clades())
    if len(clades) < 3:
        problems.append("Newick tree: fewer than 3 nodes")
        return None
    for clade in clades:
        name = clade.name or "<unnamed internal>"
        if clade.is_terminal():
            if not clade.name:
                problems.append("Newick tree: unnamed leaf")
        elif len(clade.clades) not in (2, 3):
            problems.append("Newick tree node '" + name + "': expected binary "
                            "(unrooted) branching, found " + str(len(clade.clades)) +
                            " children")
        if clade is tree.root:
            continue
        bl = clade.branch_length
        if bl is None:
            problems.append("Newick tree node '" + name + "': missing branch length")
        elif not (bl > 0 and bl < float("inf")):
            problems.append("Newick tree node '" + name + "': branch length " +
                            repr(bl) + " is not a positive finite number")
    return tree


def cross_validate(ref: Alignment, qry: Alignment, tree, problems: list[str]) -> None:
    if qry.n_cols != ref.n_cols:
        problems.append("query alignment has " + str(qry.n_cols) +
                        " columns but reference has " + str(ref.n_cols))
    for sid in sorted(set(ref.ids) & set(qry.ids)):
        problems.append("id '" + sid + "' appears in both reference and query sets")
    leaf_names = [c.name for c in tree.get_terminals()]
    if len(leaf_names) != len(set(leaf_names)):
        problems.append("Newick tree: duplicate leaf names")
    tree_set, ref_set = set(leaf_names), set(ref.ids)
    for name in sorted(tree_set - ref_set):
        problems.append("tree leaf '" + name + "' has no matching reference sequence")
    for sid in sorted(ref_set - tree_set):
        problems.append("reference id '" + sid + "' has no matching tree leaf")
