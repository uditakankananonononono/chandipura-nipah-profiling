#!/usr/bin/env python3
"""G3: structure mapping of CHPV variable/functional sites onto VSV templates.
Templates: 6U1X (VSV L+P, 3.0 A), 2GIC (VSV N-RNA). CHPV-VSV identity: L 61.3%, N ~51%.
Only regions with >=40% local identity used for structural claims (locked gate G3)."""
import json, os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
from Bio.PDB import MMCIFParser, PPBuilder
from Bio.Align import PairwiseAligner, substitution_matrices
from chpv_profile import ROOT, REF, load_proteins

prot = load_proteins()
ref_L, ref_N, ref_P = prot['L'][REF], prot['N'][REF], prot['P'][REF]
g2 = json.load(open(os.path.join(ROOT, 'results/g2_spectrum.json')))

p = MMCIFParser(QUIET=True)
ppb = PPBuilder()
def chain_seq(pdbid, chain_id):
    s = p.get_structure(pdbid, f'{ROOT}/data/structures/{pdbid}.cif')
    for model in s:
        for ch in model:
            if ch.id == chain_id:
                return ''.join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch)), ch
    return None, None

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10; aligner.extend_gap_score = -0.5

def map_positions(template_seq, query_seq):
    """query pos (1-based) -> template index (0-based into template_seq) or None."""
    a = aligner.align(template_seq, query_seq)[0]
    t, q = a[0], a[1]
    tpos = qpos = 0
    mapping = {}
    for tc, qc in zip(t, q):
        if qc != '-':
            qpos += 1
            if tc != '-':
                mapping[qpos] = tpos
        if tc != '-':
            tpos += 1
    return mapping, a

def window_identity(aln, q_lo, q_hi):
    """identity within query positions q_lo..q_hi (1-based)"""
    t, q = aln[0], aln[1]
    qpos = 0; m = n = 0
    for tc, qc in zip(t, q):
        if qc != '-':
            qpos += 1
            if q_lo <= qpos <= q_hi and tc != '-':
                n += 1
                if tc == qc: m += 1
    return (m / n * 100) if n else None

# ---- L mapping ----
tL_seq, tL_chain = chain_seq('6U1X', 'A')
map_L, aln_L = map_positions(tL_seq, ref_L)
# template CA coords: index into chain residues (only residues with CA)
residues = [r for r in tL_chain.get_residues() if r.id[0] == ' ']
# build mapping from sequence index to residue via aligned order
# PPBuilder sequence order == order of standard residues in chain
assert len(residues) >= len(tL_seq) - 5, (len(residues), len(tL_seq))

# locate VSV catalytic landmarks in template sequence
gdn_t = tL_seq.find('GDN')
hr_h_t = None
for i in range(1150, 1300):
    if tL_seq[i:i+2] == 'HR':
        hr_h_t = i; break
print(f'VSV template: GDN at seq idx {gdn_t} (res ~{gdn_t+1}), HR His at idx {hr_h_t} (res ~{hr_h_t+1 if hr_h_t else None})')

# CHPV landmark -> template position check (consistency)
g1 = json.load(open(os.path.join(ROOT, 'results/g1_positive_control.json')))
chpv_gdn = g1['gdn_pos']  # 703
print('CHPV GDN', chpv_gdn, '-> template idx', map_L.get(chpv_gdn), '(template GDN idx', gdn_t, ')')
print('CHPV HR His 1217 -> template idx', map_L.get(1217), '(template HR idx', hr_h_t, ')')

def ca_of(idx):
    if idx is None or idx >= len(residues): return None
    r = residues[idx]
    return r['CA'].get_vector() if 'CA' in r else None

gdn_ca, hr_ca = ca_of(gdn_t), ca_of(hr_h_t)

# per-domain CHPV-VSV identity (G3 quality metrics)
doms = [('RdRp',588,774),('Capping',856,1324),('Connector',1348,1547),('MTase',1629,1826)]
print('\nL per-domain CHPV-VSV identity:')
dom_id = {}
for nm, lo, hi in doms:
    dom_id[nm] = window_identity(aln_L, lo, hi)
    print(f'  {nm}: {dom_id[nm]:.1f}%' if dom_id[nm] else f'  {nm}: n/a')

# map all L variable sites + tdCE sites; distance to catalytic centers
tdce = [819, 978, 1658]
rows = []
for site in g2['spectrum']['L']:
    if site['n_variant'] == 0: continue
    pos = site['pos']
    tidx = map_L.get(pos)
    ca = ca_of(tidx)
    d_gdn = (ca - gdn_ca).norm() if (ca is not None and gdn_ca is not None) else None
    d_hr = (ca - hr_ca).norm() if (ca is not None and hr_ca is not None) else None
    rows.append(dict(pos=pos, entropy=site['entropy'], n_variant=site['n_variant'],
                     mapped=tidx is not None, d_to_GDN=round(d_gdn,1) if d_gdn else None,
                     d_to_HR_His=round(d_hr,1) if d_hr else None))
n_mapped = sum(1 for r in rows if r['mapped'])
print(f'\nL variable sites: {len(rows)}, mapped to template: {n_mapped}')
print('tdCE causal sites (structure):')
for s in tdce:
    tidx = map_L.get(s); ca = ca_of(tidx)
    d_gdn = (ca - gdn_ca).norm() if (ca is not None and gdn_ca is not None) else None
    d_hr = (ca - hr_ca).norm() if (ca is not None and hr_ca is not None) else None
    dom = next((nm for nm, lo, hi in doms if lo <= s <= hi), 'inter-domain')
    print(f'  L{s} ({dom}): mapped={tidx is not None} d_to_GDN={d_gdn and round(d_gdn,1)}A d_to_HR={d_hr and round(d_hr,1)}A')

# ---- N mapping ----
tN_seq, tN_chain = chain_seq('2GIC', 'A')
map_N, aln_N = map_positions(tN_seq, ref_N)
residues_N = [r for r in tN_chain.get_residues() if r.id[0] == ' ']
# RNA chain: find hetero/residues not standard in other chains; 2GIC RNA is chain R or similar
s2 = p.get_structure('2GIC', f'{ROOT}/data/structures/2GIC.cif')
rna_atoms = []
for model in s2:
    for ch in model:
        for r in ch:
            if r.id[0] != ' ' or ch.id in 'ABCDE':
                pass
# simpler: RNA = residues with names like A/U/G/C in any chain
for model in s2:
    for ch in model:
        for r in ch:
            if r.resname.strip() in ('A','U','G','C','U5'):
                rna_atoms.extend(list(r.get_atoms()))
print(f'\n2GIC RNA atoms: {len(rna_atoms)}')
nN_rows = []
for site in g2['spectrum']['N']:
    if site['n_variant'] == 0: continue
    pos = site['pos']
    tidx = map_N.get(pos)
    ca = None
    if tidx is not None and tidx < len(residues_N) and 'CA' in residues_N[tidx]:
        ca = residues_N[tidx]['CA'].get_vector()
    d_rna = min(( (ca - a.get_vector()).norm() for a in rna_atoms[:400]), default=None) if ca is not None else None
    nN_rows.append(dict(pos=pos, entropy=site['entropy'], n_variant=site['n_variant'],
                        mapped=tidx is not None, d_to_RNA=round(d_rna,1) if d_rna else None))
print(f'N variable sites: {len(nN_rows)}, mapped: {sum(1 for r in nN_rows if r["mapped"])}')

json.dump(dict(L_domain_identity=dom_id, L_var_sites=rows, N_var_sites=nN_rows,
               template_L_gdn_idx=gdn_t, template_L_hr_idx=hr_h_t, tdce=tdce),
         open(os.path.join(ROOT, 'results/g3_structure.json'), 'w'))
print('wrote results/g3_structure.json')
