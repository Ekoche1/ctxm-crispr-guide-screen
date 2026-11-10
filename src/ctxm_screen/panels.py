"""Build target and non-target sequence panels from the AMRFinderPlus database.

Source: NCBI Pathogen Detection "Reference Gene Catalog" and ``AMR_CDS.fa``.
Each FASTA header looks like::

    WP_000239590.1|NG_048935.1|1|1|blaCTX-M-15|blaCTX-M|product NG_048935.1:101-976

so field 4 is the allele name and field 1 the RefSeq nucleotide accession.

Design decisions (explained in docs/METHODS.md):

* The catalog has no CTX-M "group" label, so groups are *derived* from the
  data: each allele is aligned to one anchor allele per historical group and
  assigned to its best match.
* An allele is assigned only if it is a clear match (identity >= MIN_IDENTITY)
  and clearly better than the runner-up (lead >= MIN_LEAD). Hybrids and
  odd alleles go to "unassigned" and are reported, not forced into a group.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import Align

from .sequences import Record, validate_records

#: Anchor allele for each historical CTX-M group (names from the literature;
#: the grouping itself is re-measured from sequence identity in every run).
GROUP_ANCHORS = {
    "G1": "blaCTX-M-1",
    "G2": "blaCTX-M-2",
    "G8": "blaCTX-M-8",
    "G9": "blaCTX-M-9",
    "G25": "blaCTX-M-25",
}
MIN_IDENTITY = 0.93
MIN_LEAD = 0.05

#: Other beta-lactamase / resistance families used as non-targets, and the
#: maximum number of alleles sampled from each (evenly spaced, deterministic).
OTHER_FAMILIES = ["blaTEM", "blaSHV", "blaOXA", "blaCMY", "blaEC", "blaKPC", "blaNDM"]
PER_FAMILY_CAP = 25


@dataclass(frozen=True)
class Assignment:
    allele: str
    accession: str
    group: str  # "G1", "G2", ... or "unassigned"
    best_identity: float
    runner_up_identity: float
    note: str


def _aligner() -> Align.PairwiseAligner:
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 1, -1
    al.open_gap_score, al.extend_gap_score = -3, -1
    return al


def identity(aligner: Align.PairwiseAligner, a: str, b: str) -> float:
    """Matching positions divided by the longer sequence length."""
    aln = aligner.align(a, b)[0]
    matches = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        matches += sum(x == y for x, y in zip(a[s1:e1], b[s2:e2]))
    return matches / max(len(a), len(b))


def allele_name(rec: Record) -> str:
    return rec.id.split("|")[4]


def accession(rec: Record) -> str:
    return rec.id.split("|")[1]


def index_by_allele(records: list[Record]) -> dict[str, Record]:
    return {allele_name(r): r for r in records}


def families_from_catalog(catalog_path: str | Path) -> dict[str, list[str]]:
    """Map gene_family -> sorted allele names, from ReferenceGeneCatalog.txt."""
    fam: dict[str, list[str]] = {}
    with open(catalog_path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            fam.setdefault(row["gene_family"], []).append(row["allele"])
    return {k: sorted(v) for k, v in fam.items()}


def assign_groups(ctxm: list[Record], anchors: dict[str, Record]) -> list[Assignment]:
    """Assign each CTX-M allele to its closest anchor group, or 'unassigned'."""
    al = _aligner()
    out: list[Assignment] = []
    for rec in ctxm:
        scores = {g: identity(al, rec.seq, a.seq) for g, a in anchors.items()}
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        (g1, s1), (_, s2) = ranked[0], ranked[1]
        if s1 < MIN_IDENTITY:
            group, note = "unassigned", f"best identity {s1:.3f} < {MIN_IDENTITY}"
        elif s1 - s2 < MIN_LEAD:
            group, note = "unassigned", f"lead over runner-up {s1 - s2:.3f} < {MIN_LEAD}"
        else:
            group, note = g1, ""
        out.append(Assignment(allele_name(rec), accession(rec), group, s1, s2, note))
    return out


def sample_evenly(items: list[str], cap: int) -> list[str]:
    """Pick up to ``cap`` items evenly spaced through a sorted list."""
    items = sorted(items)
    if len(items) <= cap:
        return items
    step = len(items) / cap
    return [items[int(i * step)] for i in range(cap)]


def build_panels(
    cds_records: list[Record], catalog_path: str | Path
) -> dict[str, object]:
    """Return panels, assignments and exclusions. Pure function; no file I/O."""
    by_allele = index_by_allele(cds_records)
    fams = families_from_catalog(catalog_path)

    ctxm_names = [a for a in fams["blaCTX-M"] if a in by_allele]
    ctxm_valid, ctxm_excluded = validate_records([by_allele[a] for a in ctxm_names])
    anchors = {g: by_allele[a] for g, a in GROUP_ANCHORS.items()}
    assignments = assign_groups(ctxm_valid, anchors)
    group_of = {a.allele: a.group for a in assignments}

    target = [r for r in ctxm_valid if group_of[allele_name(r)] == "G1"]
    other_ctxm = [
        r for r in ctxm_valid if group_of[allele_name(r)] not in ("G1", "unassigned")
    ]
    unassigned = [r for r in ctxm_valid if group_of[allele_name(r)] == "unassigned"]

    other_bla: list[Record] = []
    other_excluded = []
    for fam in OTHER_FAMILIES:
        chosen = [by_allele[a] for a in sample_evenly(
            [a for a in fams.get(fam, []) if a in by_allele], PER_FAMILY_CAP)]
        ok, bad = validate_records(chosen)
        other_bla.extend(ok)
        other_excluded.extend(bad)

    return {
        "target_g1": target,
        "nontarget_other_ctxm_groups": other_ctxm,
        "nontarget_other_beta_lactamases": other_bla,
        "unassigned_ctxm": unassigned,
        "assignments": assignments,
        "exclusions": ctxm_excluded + other_excluded,
    }
