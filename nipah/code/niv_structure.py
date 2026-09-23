#!/usr/bin/env python3
"""G3: map G2 variable sites onto real NiV structures.
Structures: 9GJU (RNA-bound L-P, 2.8 A, L chain A res 6-2244, P chains B-E
477-709), 4CO6 (N core + P peptide), 4N5B/6EB9 (P multimerization domains).
Domain map (PMC11615333, PMC11885841):
  L: RdRp 1-969 | PRNTase 970-1452 | CD 1453-1758 | MTase 1759-2080 |
     CTD 2081-2244 (MTase/CTD boundary approximate - documented caveat)
  P: NTD 1-469 (disordered; 1-50 in 4CO6) | OD 470-578 | XD 652-709
  N: core 1-405 (4CO6 21-386) | tail 406-532 (disordered)
RSA via ShrakeRupley; core < 0.25 <= surface. Distances to catalytic sites.
"""
import json, os, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser, ShrakeRupley

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RES = os.path.join(HERE, '..', 'results')

L_DOMAINS = [(1, 969, 'RdRp'), (970, 1452, 'PRNTase'), (1453, 1758, 'CD'),
             (1759, 2080, 'MTase'), (2081, 2244, 'CTD')]
P_DOMAINS = [(1, 469, 'NTD_disordered'), (470, 578, 'OD'), (579, 651, 'OD-XD_linker'),
             (652, 709, 'XD')]
N_DOMAINS = [(1, 405, 'core'), (406, 532, 'tail_disordered')]

def domain_of(domains, site):
    for lo, hi, name in domains:
        if lo <= site <= hi:
            return name
    return 'unassigned'

def rsa_map(cif, chain_id):
    p = MMCIFParser(QUIET=True)
    s = p.get_structure('x', os.path.join(DATA, 'structures', cif))
    model = s[0]
    sr = ShrakeRupley(n_points=100)
    sr.compute(model, level='R')
    out = {}
    for r in model[chain_id]:
        if r.id[0] == ' ' and 'CA' in r:
            out[r.id[1]] = {'rsa': round(r.sasa, 4), 'ca': list(r['CA'].coord)}
    return out, s

def max_rsa(rsa, seqlen):
    """normalize by ShrakeRupley max acc per residue type is complex; use
    fixed empirical max CA-accessible SASA ~ 210 A^2 (documented choice)"""
    return {k: min(v['rsa'] / 210.0, 1.0) for k, v in rsa.items()}

spec = json.load(open(os.path.join(RES, 'g2_spectrum.json')))
out = {'structures': {'L': '9GJU chain A', 'P': '9GJU chain B + 4N5B/6EB9',
                      'N': '4CO6'},
       'domain_maps': {'L': L_DOMAINS, 'P': P_DOMAINS, 'N': N_DOMAINS},
       'rsa_normalization': 'SASA/210 fixed max (documented)',
       'proteins': {}}

print('computing SASA 9GJU (L+P complex)...')
rsa_L, s9 = rsa_map('9GJU.cif', 'A')
rsa_P, _ = rsa_map('9GJU.cif', 'B')
print('L modeled:', len(rsa_L), 'P modeled:', len(rsa_P))
print('computing SASA 4CO6 (N)...')
rsa_N, s4 = rsa_map('4CO6.cif', 'A')
print('N modeled:', len(rsa_N))

nrsa = {'L': max_rsa(rsa_L, 2244), 'P': max_rsa(rsa_P, 709), 'N': max_rsa(rsa_N, 532)}
model = s9[0]
cat_atoms = []
for pos in [831, 832, 833, 834]:
    if pos in rsa_L:
        cat_atoms.append(rsa_L[pos]['ca'])
cat_L = {'GDNE': cat_atoms}

import math
def dist(a, b):
    return math.dist(a, b)

stats = {}
for g in 'LPN':
    rows = []
    domains = {'L': L_DOMAINS, 'P': P_DOMAINS, 'N': N_DOMAINS}[g]
    rsa = {'L': rsa_L, 'P': rsa_P, 'N': rsa_N}[g]
    nr = nrsa[g]
    for s_ in spec['proteins'][g]['sites']:
        site = s_['site']
        rec = {'site': site, 'ref': s_['ref'], 'n_carriers': s_['n_carriers'],
               'host_groups': s_['host_groups'],
               'domain': domain_of(domains, site)}
        if site in nr:
            rec['rsa'] = round(nr[site], 3)
            rec['location'] = 'surface' if nr[site] >= 0.25 else 'core'
            if g == 'L' and site in rsa_L and cat_L['GDNE']:
                rec['dist_to_GDNE_A'] = round(min(dist(rsa_L[site]['ca'], c) for c in cat_L['GDNE']), 1)
        else:
            rec['location'] = 'unmodeled'
        rows.append(rec)
    out['proteins'][g] = rows
    loc = {}
    dom = {}
    for r in rows:
        loc[r['location']] = loc.get(r['location'], 0) + 1
        dom[r['domain']] = dom.get(r['domain'], 0) + 1
    stats[g] = {'by_location': loc, 'by_domain': dom}

out['stats'] = stats
json.dump(out, open(os.path.join(RES, 'g3_structure.json'), 'w'), indent=1)
print('=== G3 stats ===')
for g in 'LPN':
    print(g, 'location:', stats[g]['by_location'], 'domains:', stats[g]['by_domain'])
