#!/usr/bin/env python3
"""G5: distance from every L residue CA to bound template/product RNA in 9GJU
(chains F/G); variable vs invariant comparison. Also P/N structure panels data.
"""
import json, os, math, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
from scipy.stats import mannwhitneyu
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
struct = json.load(open(os.path.join(RES, 'g3_structure.json')))
p = MMCIFParser(QUIET=True)
s = p.get_structure('x', os.path.join(HERE, '..', 'data', 'structures', '9GJU.cif'))
m = s[0]
rna_pts = [a.coord for ch in ('F','G') for r in m[ch] for a in r if a.name in ("P","C4'","C3'")]
ca_L = {r.id[1]: r['CA'].coord for r in m['A'] if r.id[0]==' ' and 'CA' in r}
def mind(pt):
    return min(math.dist(pt, q) for q in rna_pts)
dist = {pos: round(mind(c),1) for pos, c in ca_L.items()}
var_sites = {r['site'] for r in struct['proteins']['L']}
dv = [d for pos,d in dist.items() if pos in var_sites]
di = [d for pos,d in dist.items() if pos not in var_sites]
U, pv = mannwhitneyu(dv, di, alternative='greater')
out = {'rna_chains': ['F','G'], 'n_rna_ref_points': len(rna_pts),
       'median_dist_variable': float(np.median(dv)), 'median_dist_invariant': float(np.median(di)),
       'mannwhitney_p_variable_farther': float(pv), 'n_var': len(dv), 'n_inv': len(di)}
json.dump(out, open(os.path.join(RES,'g5_rna_distance.json'),'w'), indent=1)
print(out)
# figure data dump
json.dump({str(k): v for k,v in dist.items()}, open(os.path.join(RES,'g5_L_rna_dist_per_res.json'),'w'))
