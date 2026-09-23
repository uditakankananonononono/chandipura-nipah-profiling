# chandipura-nipah-profiling

Structural + mutation profiling of Chandipura and Nipah virus replication proteins.

## Slices
- `chandipura/` - builder 4: Chandipura virus L/P/N structural + mutation profiling
- `nipah/` - builder 5: Nipah virus L/P/N structural + mutation profiling (complete)

## chandipura slice
- `chandipura/LOCKED_GATES.md` - success gates locked before outcome data
- `chandipura/code/` - chpv_profile.py (G1 positive control), chpv_spectrum.py (G2), chpv_structure.py (G3)
- `chandipura/data/` - manifests + checksums (bulk data NOT committed; regenerate via code)
- `chandipura/results/` - G1/G2/G3 JSON outputs, figures, tables
- `chandipura/paper/` - research paper (PDF + source)

## nipah slice
- `nipah/LOCKED_GATES.md` - success gates locked 2026-09-23 before outcome data
- `nipah/code/` - fetch_genomes.py (G0), niv_profile.py (G1), niv_spectrum.py (G2), niv_structure.py (G3), niv_stats.py (G4), niv_rna_distance.py (G5), niv_drug_interface.py (G6), niv_figures.py
- `nipah/data/` - manifest.json (136 genomes, per-accession md5; bulk regenerable via code)
- `nipah/results/` - g1-g6 JSON outputs, figures/, MANIFEST.sha256 (in nipah/)
- `nipah/paper/` - 21-page research paper (PDF + generator)
