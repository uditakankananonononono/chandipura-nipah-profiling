#!/usr/bin/env python3
"""chpv-profile: Chandipura virus L/P/N mutation-spectrum + structure profiler.
Reference coordinate system: strain I653514 (KF468775.1) - matches UniProt
P13179 (L) exactly; N differs from P11211 only at R37K (documented)."""
import json, os, math, hashlib
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = 'KF468775.1'  # I653514
HUMAN = ['GU190711.1', 'GU212856.1', 'GU212857.1', 'GU212858.1', 'PQ185534.2']
VECTOR = [f'MT0196{i:02d}.1' for i in range(8, 20)] + ['ON158116.1', 'ON158117.1', 'ON158118.1', 'ON158119.1', 'HM627187.1']
HEDGEHOG = ['HM627186.1']
LAB = ['KF468772.1', 'KF468773.1', 'KF468774.1']  # CH256/CH157/CH112 tdCE mutants
EXCLUDED_DUP = ['NC_020805.1']  # RefSeq identical copy of GU212856.1

def load_proteins():
    return json.load(open(os.path.join(ROOT, 'data/proteins/proteins.json')))

def group_of(acc):
    if acc in HUMAN: return 'human'
    if acc in VECTOR: return 'vector'
    if acc in HEDGEHOG: return 'hedgehog'
    if acc in LAB: return 'lab-mutant'
    if acc == REF: return 'reference-lab'
    return 'unknown'

# ---------------- curated landmarks (I653514 numbering, 1-based) -------------
LANDMARKS = {
 'L': [
   dict(name='RdRp catalytic domain', kind='region', start=588, end=774, source='UniProt P13179'),
   dict(name='RdRp motif C GDN', kind='motif', pattern='GDN', within=(588, 774), source='UniProt P13179 / NNS-virus canon'),
   dict(name='Capping domain', kind='region', start=856, end=1324, source='UniProt P13179'),
   dict(name='PRNTase domain', kind='region', start=1071, end=1321, source='UniProt P13179'),
   dict(name='priming-capping loop', kind='region', start=1142, end=1179, source='UniProt P13179'),
   dict(name='HR motif His (PRNTase nucleophile)', kind='residue', pos=1217, aa='H', tol=20, source='Ogino 2010 doi:10.1099/vir.0.019307-0; UniProt P13179 active site'),
   dict(name='HR motif Arg', kind='residue', pos=1218, aa='R', tol=20, source='Ogino 2010'),
   dict(name='Connector domain', kind='region', start=1348, end=1547, source='UniProt P13179'),
   dict(name='MTase domain', kind='region', start=1629, end=1826, source='UniProt P13179'),
   dict(name='MTase K (cat)', kind='residue', pos=1640, aa='K', tol=20, source='UniProt P13179 active site'),
   dict(name='MTase D (cat)', kind='residue', pos=1751, aa='D', tol=20, source='UniProt P13179 active site'),
   dict(name='MTase K (cat)', kind='residue', pos=1784, aa='K', tol=20, source='UniProt P13179 active site'),
   dict(name='MTase E (cat)', kind='residue', pos=1822, aa='E', tol=20, source='UniProt P13179 active site'),
   dict(name='L-P interaction 694', kind='residue', pos=694, aa=None, tol=20, source='UniProt P13179'),
   dict(name='L-P interaction 1409', kind='residue', pos=1409, aa=None, tol=20, source='UniProt P13179'),
   dict(name='L-P interaction 2007', kind='residue', pos=2007, aa=None, tol=20, source='UniProt P13179'),
   dict(name='L-P interaction 2080', kind='residue', pos=2080, aa=None, tol=20, source='UniProt P13179'),
 ],
 'P': [
   dict(name='W135 dimerization', kind='residue', pos=135, aa='W', tol=20, source='RSC Adv 2015 doi:10.1039/c5ra20863g'),
   dict(name='disordered region 1', kind='region', start=24, end=47, source='UniProt E3T2G5'),
   dict(name='acidic region', kind='region', start=55, end=69, source='UniProt E3T2G5'),
   dict(name='disordered region 3', kind='region', start=171, end=209, source='UniProt E3T2G5'),
 ],
 'N': [
   dict(name='N0-P interaction region', kind='region', start=1, end=180, source='Mondal 2012 doi:10.1371/journal.pone.0034623'),
   dict(name='N-RNA-P interaction region', kind='region', start=320, end=390, source='Mondal 2012'),
   dict(name='oligomerization region', kind='region', start=180, end=264, source='Mondal 2012'),
   dict(name='disorder 117-125', kind='region', start=117, end=125, source='Sci Rep 2021 doi:10.1038/s41598-021-92581-6'),
   dict(name='disorder 355-371', kind='region', start=355, end=371, source='Sci Rep 2021'),
 ],
}

DOMAINS_L = [  # for domain-level statistics
 ('RdRp', 588, 774), ('Capping', 856, 1324), ('Connector', 1348, 1547), ('MTase', 1629, 1826)]

def entropy(counts):
    n = sum(counts.values())
    if n == 0: return 0.0
    h = 0.0
    for c in counts.values():
        if c: p = c / n; h -= p * math.log2(p)
    return h

def build_alignment(gene, prot, accessions):
    """Return dict acc -> aligned seq (reference length, gap='-'). Only P needs realignment (HM627186 310 aa)."""
    ref = prot[gene][REF]
    out = {}
    for acc in accessions:
        s = prot[gene][acc]
        if len(s) == len(ref):
            out[acc] = s
        else:
            from Bio.Align import PairwiseAligner, substitution_matrices
            al = PairwiseAligner()
            al.mode = 'global'
            al.substitution_matrix = substitution_matrices.load('BLOSUM62')
            al.open_gap_score = -10; al.extend_gap_score = -0.5
            a = al.align(ref, s)[0]
            # project: for each ref position, aligned residue or '-'
            mapping = {}
            rpos = spos = 0
            ra, sa = a[0], a[1]
            proj = []
            for rc, sc in zip(ra, sa):
                if rc != '-':
                    proj.append(sc if sc != '-' else '-')
            out[acc] = ''.join(proj)
    return out

def run_g1(prot):
    """Positive control: landmark recovery on reference + conservation scan."""
    results = []
    for gene, lms in LANDMARKS.items():
        ref = prot[gene][REF]
        for lm in lms:
            rec = dict(gene=gene, name=lm['name'], kind=lm['kind'], source=lm['source'])
            if lm['kind'] == 'residue':
                pos, aa = lm['pos'], lm['aa']
                actual = ref[pos-1]
                rec['expected_pos'] = pos; rec['expected_aa'] = aa; rec['actual_aa'] = actual
                rec['recovered'] = (aa is None) or (actual == aa)
            elif lm['kind'] == 'motif':
                lo, hi = lm['within']
                window = ref[lo-1:hi]
                idx = window.find(lm['pattern'])
                rec['found_pos'] = (lo + idx) if idx >= 0 else None
                rec['recovered'] = idx >= 0
            else:
                rec['start'], rec['end'], rec['recovered'] = lm['start'], lm['end'], True
            results.append(rec)
    n_res = [r for r in results if r['kind'] in ('residue', 'motif')]
    n_pass = sum(1 for r in n_res if r['recovered'])
    gate = n_pass / len(n_res) >= 0.8 if n_res else False
    return results, n_pass, len(n_res), gate

def invariance_control(prot):
    """Known-answer control: catalytic residues must be invariant across all 28 genomes."""
    checks = [('L', 1217), ('L', 1218), ('L', 1640), ('L', 1751), ('L', 1784), ('L', 1822), ('P', 135)]
    refseq = {g: prot[g][REF] for g in ['L', 'P', 'N']}
    # locate GDN in reference L
    lo, hi = 588, 774
    idx = refseq['L'][lo-1:hi].find('GDN')
    gdn_pos = lo + idx
    for i in range(3):
        checks.append(('L', gdn_pos + i))
    all_accs = [a for a in prot['L'] if a not in EXCLUDED_DUP]
    aln = {g: build_alignment(g, prot, all_accs) for g in ['L', 'P', 'N']}
    rows = []
    for gene, pos in checks:
        refaa = refseq[gene][pos-1]
        variants = defaultdict(list)
        for acc in all_accs:
            aa = aln[gene][acc][pos-1]
            if aa != refaa:
                variants[aa].append(acc)
        rows.append(dict(gene=gene, pos=pos, ref_aa=refaa, invariant=len(variants) == 0,
                         variants={k: v for k, v in variants.items()}))
    return gdn_pos, rows

if __name__ == '__main__':
    prot = load_proteins()
    results, n_pass, n_tot, gate = run_g1(prot)
    print('=== G1 POSITIVE CONTROL: landmark recovery (reference I653514) ===')
    for r in results:
        if r['kind'] == 'residue':
            flag = 'OK ' if r['recovered'] else 'FAIL'
            print(f"[{flag}] {r['gene']} {r['name']}: pos {r['expected_pos']} expected {r['expected_aa'] or 'any'} actual {r['actual_aa']}")
        elif r['kind'] == 'motif':
            flag = 'OK ' if r['recovered'] else 'FAIL'
            print(f"[{flag}] {r['gene']} {r['name']}: found at {r['found_pos']}")
        else:
            print(f"[OK ] {r['gene']} {r['name']}: {r['start']}-{r['end']} registered")
    print(f'landmark recovery: {n_pass}/{n_tot} residue/motif checks -> gate {"PASS" if gate else "FAIL"}')
    gdn_pos, rows = invariance_control(prot)
    print(f'\nGDN motif at L {gdn_pos}-{gdn_pos+2} (expect within 588-774)')
    print('=== G1b KNOWN-ANSWER CONTROL: catalytic-residue invariance across 27 genomes ===')
    allpass = True
    for r in rows:
        flag = 'OK ' if r['invariant'] else 'FAIL'
        if not r['invariant']: allpass = False
        print(f"[{flag}] {r['gene']}{r['pos']}{r['ref_aa']} invariant across all genomes: {r['invariant']} {dict(r['variants']) if r['variants'] else ''}")
    print('invariance control:', 'PASS' if allpass else 'FAIL')
    json.dump(dict(landmarks=results, gdn_pos=gdn_pos, invariance=rows,
                   gate_landmarks=gate, gate_invariance=allpass),
              open(os.path.join(ROOT, 'results/g1_positive_control.json'), 'w'), indent=2)
