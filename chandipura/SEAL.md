# Chandipura slice seal - 2026-10-10 (PROVISIONAL pending auditor delta verification)

Sealed commit: da9bf8e4dd7c560a8c18f2a39699640e652ab5b1
Root manifest: MANIFEST.sha256 (sha256 of file: 1a5c9a8d6483154df843e905510d13ff5a43fdfb18c8f873a0c6b813ac3efce7, 22 entries, sha256sum -c clean)
Regeneration evidence: chandipura/G4_RUN_LOG.md (fresh-sandbox rerun of the G4 data build, 2026-10-10)

Seal criteria (all met at the sealed commit):
1. Phantom-seal wording corrected: nipah/README.md and root README.md state the slice is built and byte-verified with seal pending review; nipah/LOCKED_GATES.md stays frozen with its original text (edit 8131dd72 reverted in b28918b5) and the superseded sentence is corrected by the 2026-10-10 erratum in nipah/README.md.
2. Gate mapping reconciled: chandipura/README.md lists the actual results files (g1_positive_control, g2_spectrum, g3_structure locked; g5_temporal, g6_pairwise, g7_n_rna_distance exploratory) and documents that g5/g6/g7 map to no locked gate and are excluded from gate evidence; root README.md matches, including the nipah G5/G6 unmapped note.
3. G4 regeneration run log committed (environment, commands, outcomes, output hashes); documentary coverage 22/22 byte-checked 2026-10-10 plus logged fresh-sandbox regeneration; outcomes as measured (G1/G3/G5 byte-identical, G2 row-order-only difference, fig4 byte drift with identical inputs). No G4 PASS claim beyond what the log shows. payload_checksums.md5 corrected (code/ paths, g7, figures, proteins.json entries; regeneration-verified hashes).
4. Chronology clarification committed: gates commit 1aaa02a and results commit 2428be2 share the same author timestamp second (2026-09-23 13:40:06 +0530); the gates text's "written 13:29 IST" claim is not corroborated by the stamp. Equal author timestamps are not ordering evidence; precedence rests on DAG order (1aaa02a parents 2428be2) plus the lock attestation only.
5. Per-slice seal: this record covers the chandipura slice only. The nipah slice is NOT sealed.

Attestation artifacts: this SEAL.md and the tag seal-chandipura-2026-10-10 are attestation records created after the sealed commit. Like every manifest, the manifest cannot hash itself; SEAL.md is intentionally outside MANIFEST.sha256 and documented as such here. The tag is pushed only after the auditor delta verification passes.
