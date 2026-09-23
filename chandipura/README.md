# CHAN-STRUCT: Chandipura virus L/P/N structural + mutation profiling

Builder 4 slice of chandipura-nipah-profiling. Reference coordinate system:
strain I653514 (KF468775.1) - matches UniProt P13179 (L) exactly; N differs from
UniProt P11211 only at R37K (documented in paper).

## Layout
- `LOCKED_GATES.md` - success gates locked 2026-09-23 before outcome data
- `code/` - pipeline (Python 3.10; biopython, numpy, pandas, scipy, matplotlib)
- `data/` - manifest.json (per-accession md5), payload_checksums.md5. Bulk
  genomes/structures are NOT committed; regenerate with the fetch steps below.
- `results/` - g1-g6 JSON outputs, figures/
- `paper/` - research paper PDF + generator

## Reproduce (fresh sandbox)
1. Fetch genomes (NCBI eutils, txid11272, SLEN 10000:12000; 28 records ->
   chpv_genomes.gb / chpv_genomes_raw.fasta), parse per-accession FASTAs ->
   data/genomes/manifest.json (md5 per file). Verify against
   data/payload_checksums.md5.
2. Fetch templates: PDB 6U1X (VSV L+P, 3.0 A), 2GIC (VSV N-RNA) as mmCIF.
3. `python3 code/chpv_profile.py`   -> G1 positive controls (12/12 landmarks,
   10/10 catalytic invariance, tdCE known-answer recovery 7/7 coding + 3/3 silent)
4. `python3 code/chpv_spectrum.py`  -> G2 per-site spectrum + provenance
5. `python3 code/chpv_structure.py` -> G3 mapping onto 6U1X / 2GIC
6. `python3 code/chpv_figures.py`   -> figures fig1-fig5 + temporal analysis
7. `python3 paper/build_paper.py`   -> paper PDF

## Dataset
28 complete/near-complete CHPV genomes (10-12 kb); NC_020805.1 excluded as an
identical RefSeq copy of GU212856.1 -> 27 analyzed: 5 human (India, 2003-2024),
17 sandfly (Senegal 1992-97, Kenya 2016-17, Nigeria 1978), 1 hedgehog (Nigeria
1966), 3 lab tdCE mutants + reference I653514 (1965).
