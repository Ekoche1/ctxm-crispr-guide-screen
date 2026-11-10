"""Step 1: build the target / non-target panels and write the data manifest.

Reads  data/raw/AMR_CDS.fa and data/raw/ReferenceGeneCatalog.txt
Writes data/panels/*.fasta, data/panels/panel_manifest.tsv,
       data/panels/exclusions.tsv, data/MANIFEST.md (source hashes)
"""

import csv
import hashlib
from pathlib import Path

from ctxm_screen.panels import accession, allele_name, build_panels
from ctxm_screen.sequences import parse_fasta

ROOT = Path(__file__).resolve().parents[1]
RAW, OUT = ROOT / "data" / "raw", ROOT / "data" / "panels"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path: Path, records) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for r in records:
            fh.write(f">{allele_name(r)} {accession(r)}\n")
            for i in range(0, len(r.seq), 70):
                fh.write(r.seq[i : i + 70] + "\n")


def main() -> None:
    cds = list(parse_fasta(RAW / "AMR_CDS.fa"))
    result = build_panels(cds, RAW / "ReferenceGeneCatalog.txt")

    role = {
        "target_g1": "target",
        "nontarget_other_ctxm_groups": "non-target",
        "nontarget_other_beta_lactamases": "non-target",
        "unassigned_ctxm": "unassigned",
    }
    with open(OUT / "panel_manifest.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel", "role", "allele", "refseq_nucleotide_accession", "length_bp"])
        for panel, r_ in role.items():
            write_fasta(OUT / f"{panel}.fasta", result[panel])
            for rec in result[panel]:
                w.writerow([panel, r_, allele_name(rec), accession(rec), len(rec.seq)])

    with open(OUT / "ctxm_group_assignments.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["allele", "accession", "group", "best_identity", "runner_up_identity", "note"])
        for a in result["assignments"]:
            w.writerow([a.allele, a.accession, a.group, f"{a.best_identity:.4f}",
                        f"{a.runner_up_identity:.4f}", a.note])

    with open(OUT / "exclusions.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["record_id", "reason"])
        for e in result["exclusions"]:
            w.writerow([e.id, e.reason])

    version = (RAW / "version.txt").read_text().strip()
    lines = [
        "# Data manifest",
        "",
        "Source: NCBI Pathogen Detection, AMRFinderPlus database "
        f"(version {version}).",
        "URL: https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/",
        "",
        "| File | SHA-256 |",
        "|---|---|",
    ]
    for name in ("AMR_CDS.fa", "ReferenceGeneCatalog.txt", "version.txt"):
        lines.append(f"| {name} | `{sha256(RAW / name)}` |")
    lines += ["", "Panel sizes written by scripts/01_build_panels.py:", "",
              "| Panel | Role | Sequences |", "|---|---|---|"]
    for panel, r_ in role.items():
        lines.append(f"| {panel} | {r_} | {len(result[panel])} |")
    lines.append(f"\nExcluded records: {len(result['exclusions'])} (see data/panels/exclusions.tsv)")
    (ROOT / "data" / "MANIFEST.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    for panel in role:
        print(f"{panel:36s} {len(result[panel]):4d}")
    print("exclusions:", [(e.id.split('|')[4] if '|' in e.id else e.id, e.reason) for e in result["exclusions"]])


if __name__ == "__main__":
    main()
