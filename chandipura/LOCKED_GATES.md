# CHAN-STRUCT: Chandipura virus L/P/N structural + mutation profiling
## Locked success gates (written 2026-09-23 13:29 IST, BEFORE any outcome data is inspected)
Builder 4 slice of chandipura-nipah-profiling. Scope: Chandipura virus only (builder 5 takes Nipah).

## Research question
Do CHPV replication proteins (L polymerase, P phosphoprotein, N nucleoprotein)
carry host-stratified (human vs vector/animal) amino-acid substitution patterns,
and where do variable sites fall relative to functional domains and 3D structure?

## Hypothesis
Human-derived isolates carry a distinct set of substitutions vs vector/animal
isolates, concentrated OUTSIDE catalytic/binding cores (purifying selection on
function); human-associated substitutions, if any, map to surface/peripheral
regions of L/P/N structures.

## Prior-art verdict: CROWDED (gap identified)
Closest works:
1. PINSA 2026, doi:10.1007/s43538-026-00906-8 - AlphaFold3 models of 5 CHPV
   proteins + drug docking; 29-genome per-gene average identity only. No per-site
   mutation spectrum, no host stratification, no accession provenance.
2. IJMS 2025, doi:10.3390/ijms26031021 - phylogenetics of 23 genomes
   (human/sandfly/hedgehog); clades only, no per-residue spectrum, no structure.
3. Microbiol Spectr 2025, doi:10.1128/spectrum.01578-25 - one 2024 Gujarat genome.
4. PLoS ONE 2012, doi:10.1371/journal.pone.0030315 - whole-genome comparison.
5. Ogino 2010, doi:10.1099/vir.0.019307-0 - HR motif / PRNTase biochemistry (L).
6. tdCE host-range paper, doi:10.1099/vir.0.059204-0 - L mutations set host range.
GAP: no study combines curated-domain positive-control validation, a per-site
accession-provenanced mutation spectrum across ALL complete genomes stratified
by host, and structure-mapped interpretation of that spectrum. That is the claim.

## GATES
- G0 byte-lock: all CHPV genomes txid11272 SLEN 10000-12000 nt fetched from NCBI
  eutils; manifest of accession, host, country, year, md5 per file.
  PASS: >=20 genomes total; >=5 non-human isolates (else host-stratified claims
  downgraded and reported thin).
- G1 positive control (landmark recovery): automated scan of reference L/P/N must
  recover curated knowns BEFORE any novel claim:
  L: RdRp motifs A-F incl. catalytic GDN in motif C; HR capping motif (His-Arg);
     C-terminal MTase catalytic tetrad region (K-D-K-E).
  P: dimerization region around W135; N0-P binding region; disordered N-term.
  N: N0-P interaction region (aa 1-180); N-RNA-P region (aa 320-390).
  PASS: >=80% of curated landmarks recovered within +/-20 aa of curated position.
- G2 mutation spectrum: per-site amino-acid variability for L, P, N across all
  genomes; every substitution traceable to accession(s) in a provenance table;
  human-vs-animal contingency with Fisher exact (any cell <5); verified counts
  (script totals == table totals, printed in report).
  PASS: 100% of reported substitutions have accession-level provenance; totals
  reconcile exactly.
- G3 structure mapping: variable/conserved sites mapped onto structure.
  AlphaFold DB model for N (P11211) if pLDDT>=70 in mapped regions; L (2092 aa)
  via VSV L template PDB 6U1X homology mapping if no full-length AF model meets
  quality bar; P treated as partly disordered (Sci Rep 2021 disorder analysis),
  mapped regions restricted accordingly.
  PASS: structural claims only over validated regions (pLDDT>=70 or template
  identity >=40%); quality metrics printed per mapped region.
- G4 reproducibility: figures regenerate from manifest+code; md5 manifest covers
  every payload file; README reproduces end-to-end in fresh sandbox.

## Negative-result protocol
Null findings (e.g., no host stratification) are reported as informative
negatives with power statement, never re-fished.

## Tool deliverable
Reusable CLI: chpv-profile (download -> align -> per-site spectrum -> provenance
-> structure map -> figures), packaged in repo with README.

---
Fleet-QC note (2026-09-23 13:39 IST): re-committed standalone per orchestrator
rule "gates doc in its own commit before results commits". Original lock
timestamp stands: 2026-09-23 13:29 IST, before any outcome data was inspected
(attested in milestone-1 report to parent at 13:29 IST).
