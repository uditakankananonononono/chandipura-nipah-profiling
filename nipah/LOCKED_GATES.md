# NIV-STRUCT: Nipah virus L/P/N structural + mutation profiling
## Locked success gates (written 2026-09-23 ~18:50 IST, BEFORE any outcome data is inspected)
Builder 5 slice of chandipura-nipah-profiling. Scope: Nipah virus only
(builder 4's Chandipura slice is sealed and untouched).

## Research question
Do Nipah virus replication proteins (L polymerase, P phosphoprotein,
N nucleoprotein) carry host-stratified (human vs bat reservoir vs pig
intermediate) amino-acid substitution patterns, and where do variable sites
fall relative to functional domains and 3D structure?

## Hypothesis
Bat-reservoir isolates carry a distinct substitution set vs human isolates;
human-lineage substitutions concentrate OUTSIDE catalytic/binding cores
(purifying selection on function) and map to surface/peripheral regions of
L/P/N structures. P, known to be the most variable paramyxovirus replication
protein, should show the highest per-site variability; L the most constrained.

## Prior-art verdict: CROWDED (gap identified)
Closest works:
1. Nat Commun 2024 (cryo-EM NiV L-P polymerase complex, 3.19 A) - structure
   and drug resistance only; no mutation spectrum, no host stratification.
2. IJM 2024, doi:10.18502/ijm.v16i1.14879 - comparative genomics of human NiV
   isolates, SE Asia; variant lists, no structure mapping, no bat-pig
   stratification, no accession-level per-site provenance.
3. Virus Evol 2021, doi:10.1093/ve/veaa062 - phylodynamics 1999-2015; clades
   and rates only, no per-residue spectrum, no structure.
4. bioRxiv 2023, doi:10.1101/2023.07.14.23292668 - genetic diversity across
   spatial scales; phylogenetics only.
5. Yabukarski 2014 (NSMB, PDB 4CO6) - N0-P complex crystal structure.
6. Bruhn 2014 (JVI, PDB 4N5B) + PDB 6EB9 - P tetramerization + X domains.
GAP: no study combines curated-domain positive-control validation, a per-site
accession-provenanced mutation spectrum across ALL complete NiV genomes
stratified by host, and structure-mapped interpretation of that spectrum on
the real NiV L-P/N-P structures. That is the claim.

## GATES
- G0 byte-lock: all Nipah virus genomes txid121791, SLEN 17000-19000 nt,
  fetched from NCBI nuccore eutils; manifest of accession, host, country,
  year, md5 per file; bulk FASTA regenerable, not committed.
  PASS: >=100 genomes total; >=50 human; >=20 bat (Pteropus) isolates
  (else host-stratified claims downgraded and reported thin).
  Census pre-check (design input, not outcome data): 136 records, 81 human,
  27 Pteropus, 5 Sus scrofa as of 2026-09-23.
- G1 positive control (landmark recovery): automated scan of reference
  NC_002728.1 L/P/N must recover curated knowns BEFORE any novel claim:
  L: RdRp motifs A-F incl. catalytic GDN in motif C; HR capping motif
     (His-Arg) in CR-V; MTase catalytic K-D-K-E tetrad in CR-VI;
     domain order RdRp -> PRNTase/capping -> CR-V -> MTase -> CTD.
  P: editing (V/W) site; N-terminal N0-P binding region (4CO6); central
     tetramerization coiled-coil (4N5B); C-terminal X domain (6EB9).
  N: N-terminal arm / N0-P chaperone binding region (4CO6); conserved
     RNA-binding core; disordered C-terminal tail.
  PASS: >=90% curated landmarks recovered in expected order; catalytic
  residues invariant across the dataset (checked in G2).
- G2 mutation spectrum: per-site amino-acid spectrum of L, P, N across every
  genome vs reference, with accession-level provenance for every observed
  substitution and host stratification (human / Pteropus / pig / other).
  PASS: >=95% of genomes yield all three proteins parsed; spectrum covers
  every site of each protein; every non-reference residue traceable to
  accessions.
- G3 structure mapping: variable sites mapped onto real NiV structures
  (L-P cryo-EM PDB entry; 4CO6 N-P; 4N5B P-tet; 6EB9 P-XD; AlphaFold DB
  fallback for unresolved regions) with a per-site core/surface classification.
  PASS: mapping table for L/P/N produced; every figure regenerable from
  committed code + manifest.
- G4 honest negatives: sites/regions where host stratification is absent or
  thin are reported as such; no re-fishing after outcomes are seen.

## Methodological contribution (quantified, to be benchmarked)
An open, validated pipeline (nipah_profile/spectrum/structure) delivering
per-site, accession-provenanced, host-stratified mutation spectra - versus
per-gene average-identity (IJM 2024) and clade-only phylogenetics
(Virus Evol 2021). Benchmarks: landmark-recovery rate on positive controls;
per-site resolution vs per-gene averages; provenance completeness.
