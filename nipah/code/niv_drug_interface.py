#!/usr/bin/env python3
"""G6: conservation of known drug-relevant interfaces across all genomes.
Suramin interface (PMC11615333): L E291, K542, R551, K724, N833, K893.
GHP-88309 allosteric pocket: H1165, E922 (PMC11615333).
P-L interface residues (PMC11615333): P I576, I578, G580; L L387, K385, D384;
P XD helix a1 D657, S660, D662, R669, T670, H671.
Also substitution chemistry: Grantham distance stats per protein.
"""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
spec = json.load(open(os.path.join(RES, 'g2_spectrum.json')))

INTERFACES = {
 'suramin_interface_L': [('L', p) for p in [291, 542, 551, 724, 833, 893]],
 'GHP88309_pocket_L': [('L', p) for p in [1165, 922]],
 'P_L_contact_P': [('P', p) for p in [576, 578, 580]],
 'P_L_contact_L': [('L', p) for p in [387, 385, 384]],
 'P_XD_helix_a1_P': [('P', p) for p in [657, 660, 662, 669, 670, 671]],
}
sites = {g: {s['site']: s for s in spec['proteins'][g]['sites']} for g in 'LPN'}
out = {'interfaces': {}}
for name, plist in INTERFACES.items():
    rows = []
    for g, pos in plist:
        s = sites[g].get(pos)
        rows.append({'protein': g, 'site': pos,
                     'status': 'INVARIANT' if s is None else
                               f"VARIABLE n={s['n_carriers']} alts={ {a: len(v) for a,v in s['alts'].items()} }"})
    n_inv = sum(1 for r in rows if r['status'] == 'INVARIANT')
    out['interfaces'][name] = {'residues': rows, 'invariant': n_inv, 'total': len(rows)}
    print(name, n_inv, '/', len(rows), 'invariant')
    for r in rows:
        if r['status'] != 'INVARIANT': print('   ', r)

# Grantham distances (matrix from Grantham 1974, subset via biopython? use table)
GRANTHAM = {
 ('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,('A','Q'):91,('A','E'):107,('A','G'):60,('A','H'):86,('A','I'):94,('A','L'):96,('A','K'):106,('A','M'):84,('A','F'):113,('A','P'):27,('A','S'):99,('A','T'):58,('A','W'):148,('A','Y'):112,('A','V'):64,
 ('R','N'):86,('R','D'):96,('R','C'):180,('R','Q'):43,('R','E'):54,('R','G'):125,('R','H'):29,('R','I'):97,('R','L'):102,('R','K'):26,('R','M'):91,('R','F'):97,('R','P'):103,('R','S'):110,('R','T'):71,('R','W'):101,('R','Y'):77,('R','V'):96,
 ('N','D'):23,('N','C'):139,('N','Q'):46,('N','E'):42,('N','G'):80,('N','H'):68,('N','I'):149,('N','L'):153,('N','K'):94,('N','M'):142,('N','F'):158,('N','P'):91,('N','S'):46,('N','T'):65,('N','W'):174,('N','Y'):143,('N','V'):133,
 ('D','C'):154,('D','Q'):61,('D','E'):45,('D','G'):94,('D','H'):81,('D','I'):168,('D','L'):172,('D','K'):101,('D','M'):160,('D','F'):177,('D','P'):108,('D','S'):65,('D','T'):85,('D','W'):181,('D','Y'):160,('D','V'):152,
 ('C','Q'):154,('C','E'):170,('C','G'):159,('C','H'):174,('C','I'):198,('C','L'):198,('C','K'):202,('C','M'):196,('C','F'):205,('C','P'):169,('C','S'):112,('C','T'):149,('C','W'):215,('C','Y'):194,('C','V'):192,
 ('Q','E'):29,('Q','G'):87,('Q','H'):24,('Q','I'):129,('Q','L'):113,('Q','K'):53,('Q','M'):101,('Q','F'):116,('Q','P'):76,('Q','S'):68,('Q','T'):42,('Q','W'):130,('Q','Y'):99,('Q','V'):121,
 ('E','G'):88,('E','H'):40,('E','I'):134,('E','L'):138,('E','K'):56,('E','M'):126,('E','F'):140,('E','P'):93,('E','S'):80,('E','T'):65,('E','W'):152,('E','Y'):122,('E','V'):121,
 ('G','H'):98,('G','I'):135,('G','L'):138,('G','K'):127,('G','M'):127,('G','F'):153,('G','P'):42,('G','S'):56,('G','T'):59,('G','W'):184,('G','Y'):147,('G','V'):109,
 ('H','I'):142,('H','L'):99,('H','K'):32,('H','M'):87,('H','F'):100,('H','P'):77,('H','S'):89,('H','T'):47,('H','W'):115,('H','Y'):83,('H','V'):84,
 ('I','L'):5,('I','K'):102,('I','M'):10,('I','F'):21,('I','P'):95,('I','S'):142,('I','T'):89,('I','W'):61,('I','Y'):33,('I','V'):29,
 ('L','K'):107,('L','M'):15,('L','F'):22,('L','P'):98,('L','S'):145,('L','T'):92,('L','W'):61,('L','Y'):36,('L','V'):32,
 ('K','M'):95,('K','F'):102,('K','P'):103,('K','S'):121,('K','T'):78,('K','W'):110,('K','Y'):85,('K','V'):97,
 ('M','F'):28,('M','P'):87,('M','S'):135,('M','T'):81,('M','W'):67,('M','Y'):37,('M','V'):21,
 ('F','P'):114,('F','S'):155,('F','T'):103,('F','W'):40,('F','Y'):22,('F','V'):50,
 ('P','S'):74,('P','T'):38,('P','W'):147,('P','Y'):110,('P','V'):76,
 ('S','T'):58,('S','W'):177,('S','Y'):144,('S','V'):124,
 ('T','W'):128,('T','Y'):92,('T','V'):69,
 ('W','Y'):37,('W','V'):88,('Y','V'):55,
}
def gr(a,b):
    if a==b: return 0
    return GRANTHAM.get((a,b), GRANTHAM.get((b,a)))
chem = {}
for g in 'LPN':
    vals = []
    for s in spec['proteins'][g]['sites']:
        for alt, lst in s['alts'].items():
            d = gr(s['ref'], alt)
            if d is not None:
                vals += [d]*len(lst)
    vals = np.array(vals)
    chem[g] = {'n_substitutions': int(len(vals)), 'mean_grantham': float(vals.mean()),
               'median_grantham': float(np.median(vals)),
               'frac_conservative_lt60': float((vals<60).mean()),
               'frac_radical_ge100': float((vals>=100).mean())}
    print(g, chem[g])
out['substitution_chemistry'] = chem
json.dump(out, open(os.path.join(RES,'g6_drug_interface.json'),'w'), indent=1)
