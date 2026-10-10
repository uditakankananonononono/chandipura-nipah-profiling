# Chandipura slice seal - 2026-10-10

Sealed commit: 82c244d8e38a6feb54d7bdeabf971850f7374664
Root manifest: MANIFEST.sha256 (sha256 of file: 36a0f37a5e04ac8b12f5b14c394fb2614835e273b6ec624d846cdeca8d254fc8, 22 entries, sha256sum -c clean)
Regeneration evidence: chandipura/G4_RUN_LOG.md (fresh-sandbox rerun of the G4 data build, 2026-10-10)

Seal criteria (all met at the sealed commit):
1. nipah wording corrected: nipah/README.md and nipah/LOCKED_GATES.md state gates are built and byte-verified with seal pending review; root README.md matches.
2. Gate mapping reconciled: chandipura/README.md maps the unmapped g5/g6/g7 result files to locked gates G0-G4 and documents their exclusion from gate evidence, consistent with the root README G1/G2/G3 statements.
3. G4 regeneration run log committed (environment, commands, outcomes, output hashes); claims limited to evidence actually produced.
4. Chronology clarification committed: equal author timestamps are not precedence evidence; DAG order plus the lock attestation carries precedence.
5. Per-slice seal: this record covers the chandipura slice only. The nipah slice is NOT sealed.

Attestation artifacts: this SEAL.md and the tag seal-chandipura-2026-10-10 are attestation records created after the sealed commit. Like every manifest, the manifest cannot hash itself; SEAL.md is intentionally outside MANIFEST.sha256 and documented as such here.
