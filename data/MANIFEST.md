# Data manifest

Source: NCBI Pathogen Detection, AMRFinderPlus database (version 2026-08-07.1).
URL: https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/

| File | SHA-256 |
|---|---|
| AMR_CDS.fa | `bcc57417c5a170fb4a664de6898455990fcc296dbfe214b6874d761dad63962c` |
| ReferenceGeneCatalog.txt | `fa41ade01712841e639c91163243bace1c34755b927e69e0c8678540268a7fd4` |
| version.txt | `43ad0274317a2cbf80e01626e3ae6084a02f471be2f870d57fdbe2d35761911a` |

Panel sizes written by scripts/01_build_panels.py:

| Panel | Role | Sequences |
|---|---|---|
| target_g1 | target | 133 |
| nontarget_other_ctxm_groups | non-target | 123 |
| nontarget_other_beta_lactamases | non-target | 161 |
| unassigned_ctxm | unassigned | 11 |

Excluded records: 1 (see data/panels/exclusions.tsv)
