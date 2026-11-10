"""Step 5: genome-level pilot on 8 + 8 complete E. coli genomes.

Selection (NCBI nuccore E-utilities, run 2026-10-12; see docs/METHODS.md):
  chromosome_text_match:  Escherichia coli[Organism] AND complete genome[Title]
      AND srcdb_refseq[Properties] AND CTX-M-15[All Fields]           (8 records)
  chromosome_negative:    same, but NOT CTX-M[All Fields]             (first 8 of 4122)
  plasmid_gene_annotated: see PLASMID_QUERY, 12 chosen by a fixed rule
Text/annotation searches are only a *starting point*. Presence or absence of
the gene is decided here by sequence (ctxm_screen.genome.locate_gene).

Writes: data/genomes_manifest.tsv, results/genome_pilot_gene.tsv,
        results/genome_pilot_guides.tsv, results/genome_pilot_summary.json
"""

import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

from ctxm_screen.compare import IndexedSequence
from ctxm_screen.genome import locate_gene
from ctxm_screen.sequences import gc_fraction, parse_fasta

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "genomes"
RAW.mkdir(parents=True, exist_ok=True)
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

TEXT_MATCH_CHROMOSOMES = ["NZ_CP080117.1", "NZ_CP080116.1", "NZ_CP080119.1", "NZ_CP080118.1",
                          "NZ_CP080120.1", "NZ_CP014316.1", "NZ_CP139846.1", "NZ_CP084900.1"]
NEGATIVE_CHROMOSOMES = ["NZ_CP142783.1", "NZ_CP176678.1", "NZ_CP176687.1", "NZ_CP172867.1",
                        "NZ_CP149810.1", "NZ_CP127004.1", "NZ_CP135722.1", "NZ_CP119441.1"]

# Plasmid records whose own annotation names the gene. The first (text-match)
# query above proved unreliable: 6 of 8 hits had no CTX-M feature at all.
PLASMID_QUERY = ("Escherichia coli[Organism] AND blaCTX-M-15[Gene Name] "
                 "AND plasmid[Title] AND complete sequence[Title]")
N_PLASMIDS = 12
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def select_plasmids() -> list[str]:
    """Pick N plasmid accessions evenly spaced through the accession-sorted hit
    list (no recency bias). Saved to data/genome_selection.json and reused."""
    sel_path = ROOT / "data" / "genome_selection.json"
    if sel_path.exists():
        return json.loads(sel_path.read_text())["selected_plasmids"]

    def get(endpoint: str, **params) -> dict:
        with urllib.request.urlopen(f"{EUTILS}{endpoint}?{urllib.parse.urlencode(params)}", timeout=90) as r:
            data = json.load(r)
        time.sleep(0.5)
        return data

    found = get("esearch.fcgi", db="nuccore", term=PLASMID_QUERY, retmax=500, retmode="json")["esearchresult"]
    summ = get("esummary.fcgi", db="nuccore", id=",".join(found["idlist"]), retmode="json")["result"]
    accs = sorted(summ[u]["accessionversion"] for u in summ["uids"])
    step = len(accs) / N_PLASMIDS
    chosen = [accs[int(i * step)] for i in range(N_PLASMIDS)]
    sel_path.write_text(json.dumps({
        "query": PLASMID_QUERY, "retrieved": date.today().isoformat(),
        "total_hits_in_ncbi": int(found["count"]), "hits_considered": len(accs),
        "rule": f"{N_PLASMIDS} accessions evenly spaced through the accession-sorted list of hits considered",
        "selected_plasmids": chosen}, indent=2) + "\n", encoding="utf-8")
    return chosen


def fetch(acc: str) -> Path:
    path = RAW / f"{acc}.fasta"
    if not path.exists():
        q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text"})
        with urllib.request.urlopen(f"{EFETCH}?{q}", timeout=120) as r:
            path.write_bytes(r.read())
        time.sleep(0.5)  # stay well under NCBI's rate limit
    return path


def main() -> None:
    guides = pd.read_csv(ROOT / "results" / "guide_summary.tsv", sep="\t")
    ref = next(r for r in parse_fasta(ROOT / "data" / "panels" / "target_g1.fasta") if r.id == "blaCTX-M-15")

    manifest, gene_rows, guide_rows = [], [], []
    plasmids = select_plasmids()
    print("selected plasmids:", plasmids, flush=True)
    for role, accs in (("chromosome_text_match", TEXT_MATCH_CHROMOSOMES),
                       ("plasmid_gene_annotated", plasmids),
                       ("chromosome_negative", NEGATIVE_CHROMOSOMES)):
        for acc in accs:
            path = fetch(acc)
            recs = list(parse_fasta(path))
            assert len(recs) == 1, f"{acc}: expected 1 record, got {len(recs)}"
            rec = recs[0]
            n_nonacgt = sum(1 for c in rec.seq if c not in "ACGT")
            manifest.append({
                "accession": acc, "role": role, "description": rec.description,
                "length_bp": len(rec.seq), "gc": round(gc_fraction(rec.seq), 4),
                "non_acgt_bases": n_nonacgt, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "source_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}",
                "retrieved": date.today().isoformat(),
            })
            t0 = time.time()
            hit = locate_gene(rec.seq, ref.seq)
            gene_rows.append({"accession": acc, "role": role, **hit.__dict__})
            idx = IndexedSequence(acc, rec.seq)
            for g in guides.itertuples(index=False):
                any_mm, pam_mm = idx.best_matches(g.spacer_5to3)
                guide_rows.append({"accession": acc, "role": role, "guide": g.guide,
                                   "min_mm_any_pam": any_mm, "min_mm_with_pam": pam_mm,
                                   "exact_copies_with_pam": idx.count_exact_with_pam(g.spacer_5to3)})
            print(f"{acc} {role:18s} {len(rec.seq):>8,} bp  gene_found={hit.found} "
                  f"identity={hit.identity_over_ref:.4f}  ({time.time() - t0:.0f}s)", flush=True)

    pd.DataFrame(manifest).to_csv(ROOT / "data" / "genomes_manifest.tsv", sep="\t", index=False, lineterminator="\n")
    gene = pd.DataFrame(gene_rows)
    gene.to_csv(ROOT / "results" / "genome_pilot_gene.tsv", sep="\t", index=False, lineterminator="\n")
    gdf = pd.DataFrame(guide_rows)
    gdf.to_csv(ROOT / "results" / "genome_pilot_guides.tsv", sep="\t", index=False, lineterminator="\n")

    summary = {"run_date": date.today().isoformat(), "records": len(manifest)}
    for role, grp in gene.groupby("role"):
        summary[role] = {"n": len(grp), "gene_found_by_sequence": int(grp.found.sum())}
    (ROOT / "results" / "genome_pilot_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
