# Methods and decisions

This is an **in-silico study**. Everything here is computed from public
sequences. Nothing in this repository measures diagnostic sensitivity,
clinical specificity, limit of detection or wet-lab performance.

**Scope.** This is the first computational step of the work proposed in my accepted ASLM 2025
abstract. It covers one gene. Out of scope by design: other pathogens and genes, RPA primer
design, multiplexing, and any wet-lab or clinical claim.

## 1. Research question

> Are candidate Cas12a target sequences from *bla*CTX-M-15 sufficiently
> conserved across a defined panel of CTX-M group-1 alleles, and distinguishable
> from a defined non-target panel (other CTX-M groups and other
> β-lactamase families)? Do the surviving candidates behave consistently in real
> *E. coli* genomes and plasmids?

## 2. Why *E. coli* and *bla*CTX-M (and what that choice does and does not mean)

Evidence from PubMed, searched on 2026-10-12 (DOIs are links):

- CTX-M-15 has been detected in *E. coli* from hospital patients in Nigeria
  ([Olowe et al. 2010](https://doi.org/10.1016/j.ijantimicag.2009.10.004)) and in
  Nigerian *K. pneumoniae*
  ([Soge et al. 2006](https://doi.org/10.1093/jac/dki429)).
- A review of African surveys reports the CTX-M-1 group as dominant, with
  CTX-M-15 the main variant in humans and animals
  ([Alonso et al. 2017](https://doi.org/10.1111/lam.12724)).
- In rural Ghana, all ESBL-producing Enterobacteriaceae causing bloodstream
  infection were the CTX-M-15 genotype, and ESBL rates in neonatal *K.
  pneumoniae* bloodstream infections were 90.6%, versus 10% of *E. coli*
  ([Eibach et al. 2016](https://doi.org/10.1016/j.ijmm.2016.05.006)).
- A case report describes early-onset neonatal sepsis from ESBL-producing *E. coli*
  carrying a CTX-M gene (Naselli et al. 2022, PMID 35920879; no DOI was returned).

**Consequence.** *E. coli* + CTX-M group 1 is biologically well motivated. But
CTX-M-15 is a *gene* carried by several species (the Ghana data are dominated by
*K. pneumoniae*). A CTX-M-15 assay therefore detects the resistance gene, not
*E. coli*. The study is framed as gene-level; *E. coli* is the host context for
the genome pilot. These sources do not establish the situation in Nigerian
neonates specifically; no Nigerian neonatal data were analysed.

**Alternatives considered.** *K. pneumoniae* has the same gene and arguably
higher neonatal burden in the Ghana data, but *E. coli* has far more public
genomes and matches the abstract's first pathogen. The choice is
for feasibility and continuity with the abstract, not a claim that it is the
most clinically important pairing.

## 3. Data

Source of allele sequences: NCBI Pathogen Detection, AMRFinderPlus database
**version 2026-08-07.1** (`AMR_CDS.fa`, `ReferenceGeneCatalog.txt`). File
SHA-256 hashes are in `data/MANIFEST.md`. These are curated *reference alleles*,
one sequence per named allele, **not** isolates. Every allele counts equally
regardless of how common it is clinically.

### Panels (`scripts/01_build_panels.py`, `src/ctxm_screen/panels.py`)

| Panel | Definition | n |
|---|---|---|
| Target | CTX-M alleles assigned to group 1 | 133 |
| Non-target A | CTX-M alleles assigned to groups 2, 8, 9, 25 | 123 |
| Non-target B | Up to 25 alleles each of blaTEM, blaSHV, blaOXA, blaCMY, blaEC, blaKPC, blaNDM (evenly spaced, deterministic) | 161 |
| Unassigned | CTX-M alleles that fit no group clearly (hybrids, divergent) | 11 |
| Excluded | blaCTX-M-161 (contains ambiguous base Y) | 1 |

**Group assignment is measured, not assumed.** The catalog has no group label.
Each allele is globally aligned to one anchor allele per historical group
(CTX-M-1, -2, -8, -9, -25) and assigned to the best match only if identity is
at least 0.93 **and** it leads the runner-up by at least 0.05. Hybrids sit between
groups; forcing them into one would be arbitrary. They are kept in a separate
"unassigned" file and reported, not hidden. Thresholds are my choices, not
taken from a standard.

Between-group divergence is large: CTX-M-15 is 98.7% identical to the CTX-M-1
anchor and 74–76% identical to the group 2, 8, 9 and 25 anchors.

## 4. Candidate guides (`src/ctxm_screen/guides.py`)

Cas12a needs a **5'-TTTV-3' PAM** (V = A, C or G) immediately 5' of a 20-nt
protospacer. The crRNA spacer has the protospacer's sequence. We enumerate every
PAM-valid 20-mer on **both strands** of the CTX-M-15 reference (blaCTX-M-15,
876 bp, NG_048935.1). That gives 20 candidates (2 on the + strand, 18 on the −
strand). The skew is expected: the coding strand is rich in AAA, which reads as
TTT on the opposite strand. Spacer length 20 is a common choice; Cas12a
variants use roughly 20–24 nt.

PathoGD uses TTTN. We use the narrower, more standard TTTV.

## 5. Screening (`src/ctxm_screen/compare.py`, `scripts/02_screen_guides.py`)

For every guide and every panel sequence we count substitutions in every
20-nt window on both strands and keep the minimum, reported twice: ignoring the
PAM and requiring a valid TTTV PAM.

**Criteria fixed before looking at results:**

- *Conserved*: exact match **with** valid PAM in at least 95% of the target panel.
- *Specific*: at least 3 mismatches (valid PAM) to **every** non-target sequence.

The 3-mismatch rule mirrors PathoGD's published "at least 3 mismatches" rule.
The 95% cut-off is my choice.

## 6. Verification

1. **Unit tests** (`pytest`): 30 tests for this study. Tests show the code does what it is specified to
   do. They do not show the biology is right.
2. **Independent re-computation** (`scripts/03_independent_check.py`): every
   guide x sequence minimum-mismatch value is recomputed with a different
   engine (the `regex` library's fuzzy matching with the PAM written into the
   pattern) and compared. Guide coordinates, PAMs and candidate counts are
   re-derived with plain string operations.
3. **Gene found by sequence, not by label** (`src/ctxm_screen/genome.py`):
   in the genome pilot, presence of CTX-M-15 is decided by seed matching plus
   alignment. Record annotation is used only as a starting point, because it
   proved unreliable (see the README).

## 7. Genome-level pilot (`scripts/05_genome_pilot.py`)

Complete *E. coli* sequences from NCBI nuccore (E-utilities): 8 chromosome
records returned by a CTX-M-15 text search, 12 plasmids whose annotation names
blaCTX-M-15 (chosen by a fixed, saved rule), and 8 chromosome records returned
by a "NOT CTX-M" search. Queries, hit counts and chosen accessions are in
`data/genome_selection.json`; accession, length, GC, non-ACGT count and file
SHA-256 for every sequence are in `data/genomes_manifest.tsv`.

For each sequence we (a) locate the CTX-M-15 gene by sequence and report
identity over the full reference length, and (b) scan all 20 guides across the
whole sequence on both strands.

**Definitions used in the reported numbers.** A sequence is a *verified gene-positive*
if the located gene has at least 99% identity to the CTX-M-15 reference
(`scripts/07_report_tables.py`); all others are *gene-free*. Result: 14
gene-positive (12 plasmids + 2 chromosomes) and 14 gene-free chromosomes (8 from the
"NOT CTX-M" query + 6 from the CTX-M-15 text query that proved to carry no CTX-M gene).
Because every gene-positive carries the identical allele, exact guide matches in them
confirm the guides sit in real plasmid/chromosome contexts but are **not** independent
evidence of conservation; that comes from the allele panel (section 5).

## 7b. Follow-up tables

`scripts/06_allele_misses.py` lists which group-1 alleles each passing guide fails to
match exactly, and how many unassigned (hybrid/divergent) alleles it matches.
`scripts/07_report_tables.py` builds the summary tables and rewrites the README's
results table from the result files; `tests/test_report_consistency.py` fails if the
README table, run information or independent-check result drift apart.
`scripts/run_all.ps1` runs steps 1-7 and the tests in order (Windows).

## 8. Limitations (read these before citing anything)

- Reference alleles are not isolates; there is no weighting by clinical frequency.
- Mismatch counting treats all positions alike. Real Cas12a is much less
  tolerant of mismatches in the PAM-proximal seed, and bulges are ignored.
  Passing "at least 3 mismatches" does not prove no cleavage; failing it does
  not prove cleavage.
- Only 20-nt spacers, TTTV PAM, and a single reference (CTX-M-15) were used.
- No amplification (RPA/PCR) primers were designed or evaluated, and no
  published CTX-M Cas12a assay was found or evaluated (PubMed queries on
  2026-10-12 returned no match; this is not proof that none exists).
- Specificity was not tested against the human genome, the microbiome, or a
  wide range of non-*E. coli* species.
- The genome pilot is small and not representative of Nigerian isolates or of
  neonatal infections.
