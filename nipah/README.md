# NIV-STRUCT: Nipah virus L/P/N structural + mutation profiling

Builder 5 slice of chandipura-nipah-profiling. Reference coordinate system:
NC_002728.1 (NiV Malaysia 1999). The Chandipura slice (builder 4) is sealed
and untouched.

## Layout
- `LOCKED_GATES.md` - success gates locked 2026-09-23 18:50 IST, standalone
  commit, before any outcome data
- `code/` - pipeline (Python 3.10; biopython, numpy, pandas, scipy,
  matplotlib, reportlab, Pillow)
- `data/` - manifest.json (per-accession md5, 136 genomes, txid121791).
  Bulk genomes/structures NOT committed; regenerate via the fetch steps.
- `results/` - g1-g6 JSON outputs, figures/ (fig1-fig10)
- `paper/` - build_paper.py -> niv_struct_paper.pdf (21 pages)

## Reproduce (fresh sandbox)
1. `python3 code/fetch_genomes.py` -> NCBI eutils txid121791 SLEN
   17000:19000 (136 records) -> data/genomes + data/manifest.json (md5 per
   accession).
2. Structures: download 9GJU, 9IR3, 4CO6, 4N5B, 6EB9 mmCIFs from RCSB into
   data/structures/ (https://files.rcsb.org/download/<PDB>.cif).
3. `python3 code/niv_profile.py`   -> G1 positive control (20/21 = 95.2%)
4. `python3 code/niv_spectrum.py`  -> G2 per-site spectrum + provenance
5. `python3 code/niv_structure.py` -> G3 mapping onto 9GJU/4CO6
6. `python3 code/niv_stats.py`     -> G4 statistics
7. `python3 code/niv_rna_distance.py`   -> G5 RNA-distance analysis
8. `python3 code/niv_drug_interface.py` -> G6 interface conservation + chemistry
9. `python3 code/niv_figures.py`   -> figures
10. `python3 paper/build_paper.py` -> paper PDF

## Dataset
136 complete NiV genomes; 130 parsed (95.6%) including 12 ORF-recovered;
6 excluded (5 without CDS, 1 unrecoverable L). Hosts: 81 human, 27 Pteropus,
4 pig, 1 dog, 2 bat-unsp., 21 unknown. Ambiguous X/B/Z calls tracked
separately, never counted as substitutions.
