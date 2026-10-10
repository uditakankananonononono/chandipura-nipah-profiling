# G4 reproducibility regeneration - fresh-sandbox run log
Date: 2026-10-10 (IST). Prepared by the project's documentary builder lane as part of the 2026-10-10 seal-readiness amendments.
Scope: chandipura slice (CHAN-STRUCT). Claim: only evidence actually produced below.

## Environment
- Linux 6.1.158+ x86_64; Python 3.10.12
- biopython 1.88, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, matplotlib 3.10.9, reportlab 3.6.8
- Fresh directory: slice copied from repo working tree; no bulk data carried over
  (data/genomes, data/proteins, data/structures re-created empty).

## Commands
1. Genome set: NCBI eutils esearch db=nuccore term="txid11272[Organism:exp] AND
   10000:12000[SLEN]" (28 records today, unchanged from 2026-09-23), then
   efetch rettype=gb and rettype=fasta for those 28 UIDs.
2. Templates: https://files.rcsb.org/download/6U1X.cif , /2GIC.cif.
3. Derived payloads re-created from the fetched genbank: per-accession FASTAs
   (data/genomes/per_accession/), CDS-translation protein set
   (data/proteins/{G,L,M,N,P}.fasta, all 5 genes x 28 accessions), and
   data/proteins/proteins.json (uncommitted derived payload, format documented
   in code/chpv_profile.py:load_proteins).
4. Pipeline: python3 code/chpv_profile.py; code/chpv_spectrum.py;
   code/chpv_structure.py; code/chpv_figures.py. All four exited 0.
5. Paper: python3 paper/build_paper.py (required workaround, see F5).

## Outcomes (regenerated vs committed)
| Artifact | Result |
|---|---|
| chpv_genomes_raw.fasta | byte-identical (md5 f110ba90c7376c7abe8146bcab249c98) |
| chpv_genomes.gb | md5 differs (NCBI flatfile metadata drift); sequence content verified unchanged via the byte-identical raw FASTA |
| 6U1X.cif, 2GIC.cif | byte-identical (md5 8aba9ed8..., 86ae92e7...) |
| results/g1_positive_control.json | byte-identical |
| results/g2_spectrum.json | content-identical: subs_table 4768/4768 rows set-equal, all aggregate sections (spectrum, host_hits, human/vector_specific, l_domain_stats, counts) byte-equal; row ORDER differs (reconstructed proteins.json key order) |
| results/g3_structure.json | byte-identical |
| results/g5_temporal.json | byte-identical |
| results/g6_pairwise.json, g7_n_rna_distance.json | NOT regenerated - no committed script writes them (exploratory outputs, unmapped to locked gates; see README "Gate mapping status") |
| figures fig1,fig2,fig3,fig5 | byte-identical |
| fig4_structure_distances.png | bytes differ with byte-identical inputs (suspected plotting nondeterminism); flagged thin |
| paper/CHAN_STRUCT_paper.pdf | regenerated: 18 pages, 6568 words vs committed 18 pages, 6527 words; bytes differ (generator metadata) |

## Findings
- F1: Fetch/parse scripts for derived payloads are not committed (README prose
  only). Per-accession FASTA byte format was not recoverable (24 header/wrap
  variants tried), so per-accession md5s in data/manifest.json cannot be
  re-attested; their source IS byte-verified (raw FASTA match).
- F2: data/proteins/proteins.json is required by the pipeline but is neither
  committed nor covered by payload_checksums.md5; reconstructed here from CDS
  translations of the fetched genbank.
- F3: payload_checksums.md5 lists code under stale src/ paths (actual: code/),
  and omits g7_n_rna_distance.json, figures, and proteins.json.
- F4: g6/g7 have no committed generator - consistent with their unmapped,
  excluded-from-gate-evidence status.
- F5: paper/build_paper.py references docs/LOCKED_GATES.md (stale path);
  worked around by copying LOCKED_GATES.md to docs/ in the sandbox.

## Verdict
G4 regeneration: core pipeline (G1/G2/G3, figures, paper) reproduces committed
results from freshly fetched public data, byte-identical except as itemized
(G2 row order, fig4 bytes). Exceptions F1-F3/F5 are documentation gaps in the
reproduction kit, not result drift. Recorded for the review lane.
