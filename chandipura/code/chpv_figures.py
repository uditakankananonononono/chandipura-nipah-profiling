#!/usr/bin/env python3
"""Figures + supplementary analyses for CHAN-STRUCT."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from chpv_profile import ROOT, REF, HUMAN, VECTOR, HEDGEHOG, LAB, load_proteins, build_alignment

prot = load_proteins()
g2 = json.load(open(f'{ROOT}/results/g2_spectrum.json'))
g3 = json.load(open(f'{ROOT}/results/g3_structure.json'))
FIG = f'{ROOT}/results/figures'
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'figure.dpi': 150, 'font.size': 8})

# ---------- extra analysis: human-lineage temporal evolution ----------
ALL_H = ['GU212858.1', 'GU212856.1', 'GU190711.1', 'GU212857.1']  # 2003-2007 Indian human
ALLA = HUMAN + VECTOR + HEDGEHOG + LAB + [REF]
aln = {g: build_alignment(g, prot, ALLA) for g in ['L','P','N']}
temporal = {}
for g in ['N','P','L']:
    ref = prot[g][REF]
    # consensus of 2003-07 human isolates
    cons_subs = []   # 2024 isolate vs 2003-07 consensus
    for i in range(len(ref)):
        col = Counter(aln[g][a][i] for a in ALL_H)
        cons, n = col.most_common(1)[0]
        new = aln[g]['PQ185534.2'][i]
        if new != cons and n >= 3:
            cons_subs.append(dict(gene=g, pos=i+1, old=cons, new=new, support=n))
    temporal[g] = cons_subs
    print(f'{g}: 2024-vs-2003/07-consensus substitutions: {len(cons_subs)}')
    for s in cons_subs: print(f"   {g}{s['pos']} {s['old']}->{s['new']} (support {s['support']}/4)")
json.dump(temporal, open(f'{ROOT}/results/g5_temporal.json','w'), indent=2)

# ---------- extra: N variable sites vs RNA distance ----------
dvar = [r['d_to_RNA'] for r in g3['N_var_sites'] if r['d_to_RNA'] is not None]
print(f'\nN variable sites mean distance to RNA: {np.mean(dvar):.1f}A (n={len(dvar)})')

# ---------- Fig 1: dataset overview ----------
fig, ax = plt.subplots(figsize=(5.5, 3))
groups = {'Human (India)': HUMAN, 'Sandfly (Senegal)': VECTOR[:12], 'Sandfly (Kenya)': VECTOR[12:16],
          'Sandfly (Nigeria)': ['HM627187.1'], 'Hedgehog (Nigeria)': HEDGEHOG, 'Lab tdCE mutants': LAB + [REF]}
yrs = {'Human (India)': [2003,2004,2004,2007,2024], 'Sandfly (Senegal)': [1992,1994,1995,1995,1995,1995,1997,1997,1997,1997,1997,1997],
       'Sandfly (Kenya)': [2016,2016,2017,2017], 'Sandfly (Nigeria)': [1978], 'Hedgehog (Nigeria)': [1966],
       'Lab tdCE mutants': [1965,1965,1965,1965]}
cols = {'Human (India)':'#c0392b','Sandfly (Senegal)':'#2980b9','Sandfly (Kenya)':'#27ae60',
        'Sandfly (Nigeria)':'#8e44ad','Hedgehog (Nigeria)':'#e67e22','Lab tdCE mutants':'#7f8c8d'}
import numpy as _np
labels_order = list(yrs.keys())
for k, (gname, ys) in enumerate(yrs.items()):
    ys = _np.array(ys, dtype=float)
    # jitter identical years horizontally for visibility
    seen = {}
    jy = []
    for y in ys:
        seen[y] = seen.get(y, 0) + 1
        jy.append(y + (seen[y]-1)*0.55)
    ax.scatter(jy, [gname]*len(ys), s=45, color=cols[gname], edgecolor='k', linewidth=0.4, zorder=3)
ax.set_xlim(1962, 2027)
ax.set_yticks(range(len(labels_order)))
ax.set_yticklabels([f'{g}  (n={len(yrs[g])})' for g in labels_order], fontsize=7.5)
ax.set_xlabel('Collection year'); ax.set_title(f'CHPV complete-genome dataset (n=27 after dedup)')
ax.grid(axis='x', alpha=0.3); plt.tight_layout(); plt.savefig(f'{FIG}/fig1_dataset.png'); plt.close()

# ---------- Fig 2: entropy tracks ----------
dom_cfg = {
 'N': [(1,180,'N0-P int.'), (180,264,'oligom.'), (320,390,'N-RNA-P')],
 'P': [(24,47,'disord.'), (55,74,'disord./acidic'), (171,209,'disord.')],
 'L': [(588,774,'RdRp'), (856,1324,'Capping'), (1348,1547,'Connector'), (1629,1826,"2'O-MTase")]}
fig, axes = plt.subplots(3, 1, figsize=(7, 5.5), sharex=False)
for ax, gene in zip(axes, ['N','P','L']):
    sp = g2['spectrum'][gene]
    xs = [s['pos'] for s in sp]; ys = [s['entropy'] for s in sp]
    ax.fill_between(xs, ys, color='#2c3e50', alpha=0.75, lw=0)
    for lo, hi, nm in dom_cfg[gene]:
        ax.axvspan(lo, hi, color='#e67e22', alpha=0.15)
        ax.text((lo+hi)/2, max(ys)*1.02, nm, ha='center', fontsize=6.5, color='#d35400')
    ax.set_ylabel(f'{gene}\nentropy (bits)'); ax.set_xlim(1, len(sp))
    ax.spines[['top','right']].set_visible(False)
axes[2].set_xlabel('residue position (I653514 numbering)')
fig.suptitle('Per-site Shannon entropy across 27 CHPV genomes (domains shaded)')
plt.tight_layout(); plt.savefig(f'{FIG}/fig2_entropy_tracks.png'); plt.close()

# ---------- Fig 3: L domain conservation ----------
fig, ax = plt.subplots(figsize=(5, 3))
ds = g2['l_domain_stats']
labels = [d['domain'] for d in ds] + ['outside\ndomains']
means = [d['mean_entropy'] for d in ds]
rest = [s['entropy'] for i, s in enumerate(g2['spectrum']['L'], start=1)
        if not any(d['start'] <= i <= d['end'] for d in ds)]
means.append(sum(rest)/len(rest))
bars = ax.bar(labels, means, color=['#c0392b','#e67e22','#f1c40f','#27ae60','#95a5a6'], edgecolor='k', linewidth=0.5)
for b, m in zip(bars, means): ax.text(b.get_x()+b.get_width()/2, m+0.002, f'{m:.3f}', ha='center', fontsize=7)
ax.set_ylabel('mean per-site entropy (bits)'); ax.set_title('L protein: catalytic domains are purged of variation')
ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig(f'{FIG}/fig3_domain_conservation.png'); plt.close()

# ---------- Fig 4: structure distances ----------
fig, axes = plt.subplots(1, 2, figsize=(7, 2.8))
d = [r['d_to_GDN'] for r in g3['L_var_sites'] if r['d_to_GDN']]
axes[0].hist(d, bins=25, color='#2980b9', edgecolor='k', linewidth=0.4)
for s, c in [(819,'#c0392b'), (978,'#e67e22'), (1658,'#27ae60')]:
    td = [r['d_to_GDN'] for r in g3['L_var_sites'] if r['pos']==s and r['d_to_GDN']]
    if td: axes[0].axvline(td[0], color=c, lw=1.5, label=f'tdCE L{s}')
axes[0].set_xlabel('distance to GDN catalytic motif (A)'); axes[0].set_ylabel('variable sites')
axes[0].legend(fontsize=6.5); axes[0].set_title('L variable sites vs RdRp active site')
dN = [r['d_to_RNA'] for r in g3['N_var_sites'] if r['d_to_RNA']]
axes[1].hist(dN, bins=20, color='#8e44ad', edgecolor='k', linewidth=0.4)
axes[1].set_xlabel('distance to encapsidated RNA (A)'); axes[1].set_ylabel('variable sites')
axes[1].set_title('N variable sites vs genomic RNA')
for ax in axes: ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig(f'{FIG}/fig4_structure_distances.png'); plt.close()

# ---------- Fig 5: host-stratified provenance heatmap ----------
fig, ax = plt.subplots(figsize=(7, 3.4))
sites = g2['host_hits']
sites_f = [s for s in sites if s['fisher_p'] < 0.01]
accs = HUMAN + VECTOR
mat = np.zeros((len(accs), len(sites_f)))
for j, s in enumerate(sites_f):
    g, p = s['gene'], s['pos']
    for i, a in enumerate(accs):
        mat[i, j] = 0 if aln[g][a][p-1] == s['ref'] else 1
ax.imshow(mat, aspect='auto', cmap='Greys', interpolation='nearest')
ax.set_yticks(range(len(accs)))
ax.set_yticklabels([a.split('.')[0] for a in accs], fontsize=5.5)
ax.set_xticks([])
ax.axhline(len(HUMAN)-0.5, color='red', lw=1.2)
ax.set_xlabel(f'{len(sites_f)} host-differentiated sites (Fisher p<0.01), ordered N -> P -> L\nblack = residue differs from I653514 reference; red line separates human (above) from vector isolates', fontsize=7)
ax.set_title('Human (India) vs African vector isolates: deep lineage split', fontsize=9)
plt.tight_layout(); plt.savefig(f'{FIG}/fig5_host_split.png'); plt.close()
print('\nfigures written:', sorted(os.listdir(FIG)))
