#!/usr/bin/env python3
"""Figures fig1-fig6 from results JSONs. All reproducible."""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(RES, 'figures'); os.makedirs(FIG, exist_ok=True)
spec = json.load(open(os.path.join(RES, 'g2_spectrum.json')))
struct = json.load(open(os.path.join(RES, 'g3_structure.json')))
stats = json.load(open(os.path.join(RES, 'g4_stats.json')))
man = json.load(open(os.path.join(HERE, '..', 'data', 'manifest.json')))
meta = {g['accession']: g for g in man['genomes']}

DOM = {'L': [(1,969,'RdRp','#4C72B0'),(970,1452,'PRNTase','#55A868'),(1453,1758,'CD','#C44E52'),(1759,2080,'MTase','#8172B3'),(2081,2244,'CTD','#937860')],
       'P': [(1,469,'NTD (disordered)','#4C72B0'),(470,578,'OD','#55A868'),(579,651,'linker','#C44E52'),(652,709,'XD','#8172B3')],
       'N': [(1,405,'core','#4C72B0'),(406,532,'tail (disordered)','#55A868')]}

# fig1 dataset overview
fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
hc = Counter(x['host'] for x in spec['genomes'])
axes[0].bar(hc.keys(), hc.values(), color='#4C72B0')
axes[0].set_title('Analyzed genomes by host'); axes[0].set_ylabel('n')
years = []
for x in spec['genomes']:
    import re
    m = re.search(r'(19|20)\d{2}', str(meta[x['accession']]['collection_date']))
    if m: years.append(int(m.group()))
axes[1].hist(years, bins=20, color='#55A868'); axes[1].set_title('Collection years'); axes[1].set_ylabel('n')
for i, g in enumerate('LPN'):
    idn = [x['identity'][g] for x in spec['genomes']]
    axes[2].plot(sorted(idn), np.linspace(0,1,len(idn)), label=g)
axes[2].legend(); axes[2].set_title('Pairwise identity to ref (CDF)'); axes[2].set_xlabel('identity'); axes[2].set_ylabel('fraction')
plt.tight_layout(); plt.savefig(f'{FIG}/fig1_dataset.png', dpi=200); plt.close()

# fig2 spectrum tracks
fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=False)
for ax, g in zip(axes, 'LPN'):
    for lo, hi, name, c in DOM[g]:
        ax.axvspan(lo, hi, color=c, alpha=0.12)
        ax.text((lo+hi)/2, ax.get_ylim()[1]*0.02, name, ha='center', fontsize=8)
    sites = spec['proteins'][g]['sites']
    x = [s['site'] for s in sites]; y = [s['n_carriers'] for s in sites]
    ax.scatter(x, y, s=8, c='#333333')
    hx = [s['site'] for s in struct['proteins'][g] if set(s['host_groups'])=={'human'}]
    hy = [s['n_carriers'] for s in spec['proteins'][g]['sites'] if s['site'] in hx]
    ax.scatter(hx, hy, s=14, c='#C44E52', label='human-exclusive', zorder=3)
    bx = [s['site'] for s in struct['proteins'][g] if set(s['host_groups'])=={'bat'}]
    by = [s['n_carriers'] for s in spec['proteins'][g]['sites'] if s['site'] in bx]
    ax.scatter(bx, by, s=14, c='#55A868', label='bat-exclusive', zorder=3)
    ax.set_ylabel(f'{g}: carriers/site'); ax.set_xlim(0, spec['proteins'][g]['length']+20)
    ax.legend(fontsize=7, loc='upper right')
axes[2].set_xlabel('residue (reference NC_002728.1 coordinates)')
plt.tight_layout(); plt.savefig(f'{FIG}/fig2_spectrum_tracks.png', dpi=200); plt.close()

# fig3 RSA enrichment
fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
for ax, g in zip(axes, 'LPN'):
    e = stats['rsa_enrichment'].get(g, {})
    if 'var_surface_frac' not in e: ax.set_title(f'{g}: n/a'); continue
    ax.bar(['variable','invariant'], [e['var_surface_frac'], e['inv_surface_frac']], color=['#C44E52','#4C72B0'])
    ax.set_ylim(0,1); ax.set_title(f"{g} surface fraction (Fisher p={e['surface_fisher_p']:.3g})")
plt.tight_layout(); plt.savefig(f'{FIG}/fig3_rsa.png', dpi=200); plt.close()

# fig4 distance to GDNE for L
rows = [r for r in struct['proteins']['L'] if 'dist_to_GDNE_A' in r]
fig, ax = plt.subplots(figsize=(7,4))
ax.hist([r['dist_to_GDNE_A'] for r in rows], bins=30, color='#4C72B0')
ax.set_xlabel('CA distance to GDNE catalytic motif (A)'); ax.set_ylabel('variable sites')
ax.set_title('L variable-site distances to polymerase active site (9GJU)')
plt.tight_layout(); plt.savefig(f'{FIG}/fig4_gdne_distance.png', dpi=200); plt.close()

# fig5 domain variability rates
fig, ax = plt.subplots(figsize=(8,4))
names, vals, cols = [], [], []
for g in 'LPN':
    for lo,hi,name,c in DOM[g]:
        nm = name.replace('_disordered',' (dis.)')
        r = stats['domain_enrichment'][g]['sites_per_1k']
        key = [k for k in r if k.replace('_disordered',' (dis.)')==nm or k==name]
        k = list(r.keys())[ [i for i,(lo2,hi2,n2,c2) in enumerate(DOM[g]) if n2==name][0] ] if False else list(r.keys())[[i for i,x in enumerate(DOM[g])][0]]
    # simpler: iterate stats directly
for g in 'LPN':
    for k,v in stats['domain_enrichment'][g]['sites_per_1k'].items():
        names.append(f'{g}:{k}'); vals.append(v)
ax.barh(names, vals, color='#55A868'); ax.set_xlabel('variable sites per 1000 residues')
plt.tight_layout(); plt.savefig(f'{FIG}/fig5_domain_rates.png', dpi=200); plt.close()

# fig6 temporal
import re
xs, ys = [], []
for x in spec['genomes']:
    m = re.search(r'(19|20)\d{2}', str(meta[x['accession']]['collection_date']))
    if m: xs.append(int(m.group())); ys.append(x['n_subs_L'])
fig, ax = plt.subplots(figsize=(7,4))
ax.scatter(xs, ys, s=12, alpha=0.6, c='#4C72B0')
t = stats['temporal_L']
ax.set_title(f"L substitutions vs year (Spearman rho={t['spearman_rho']:.2f}, p={t['p']:.2g}, n={t['n']})")
ax.set_xlabel('collection year'); ax.set_ylabel('n substitutions vs ref (L)')
plt.tight_layout(); plt.savefig(f'{FIG}/fig6_temporal.png', dpi=200); plt.close()
print('figures done:', os.listdir(FIG))

# ---- additional figures (v2) ----
g5 = json.load(open(os.path.join(RES,'g5_rna_distance.json')))
dpr = json.load(open(os.path.join(RES,'g5_L_rna_dist_per_res.json')))
var_L = {r['site'] for r in struct['proteins']['L']}
dv = [v for k,v in dpr.items() if int(k) in var_L]
di = [v for k,v in dpr.items() if int(k) not in var_L]
fig, ax = plt.subplots(figsize=(7,4))
ax.hist([di, dv], bins=40, label=['invariant','variable'], color=['#4C72B0','#C44E52'], alpha=0.8)
ax.set_xlabel('CA distance to bound RNA (A, 9GJU chains F/G)'); ax.set_ylabel('residues')
ax.legend(); ax.set_title(f"L residues vs RNA distance (variable farther, p={g5['mannwhitney_p_variable_farther']:.2g})")
plt.tight_layout(); plt.savefig(f'{FIG}/fig8_rna_distance.png', dpi=200); plt.close()

# fig9: N (4CO6) + P OD (9GJU B) CA scatters
import warnings; warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
p_ = MMCIFParser(QUIET=True)
s4 = p_.get_structure('x', os.path.join(HERE,'..','data','structures','4CO6.cif'))
caN = {r.id[1]: r['CA'].coord for r in s4[0]['A'] if r.id[0]==' ' and 'CA' in r}
varN = {r['site']: r['n_carriers'] for r in struct['proteins']['N']}
s9 = p_.get_structure('y', os.path.join(HERE,'..','data','structures','9GJU.cif'))
caP = {r.id[1]: r['CA'].coord for r in s9[0]['B'] if r.id[0]==' ' and 'CA' in r}
varP = {r['site']: r['n_carriers'] for r in struct['proteins']['P']}
fig = plt.figure(figsize=(12,5))
for i,(ca,var,t) in enumerate([(caN,varN,'NiV N core (4CO6)'),(caP,varP,'NiV P OD+XD (9GJU chain B)')]):
    ax = fig.add_subplot(1,2,i+1,projection='3d')
    xs=[c[0] for c in ca.values()]; ys=[c[1] for c in ca.values()]; zs=[c[2] for c in ca.values()]
    cs=['#d62728' if s_ in var else '#4C72B0' for s_ in ca]
    ss=[8+3*min(var.get(s_,0),15) if s_ in var else 4 for s_ in ca]
    ax.scatter(xs,ys,zs,c=cs,s=ss,alpha=0.6,linewidths=0)
    ax.set_title(f'{t} - red = variable', fontsize=9); ax.set_axis_off(); ax.view_init(elev=15,azim=45)
plt.tight_layout(); plt.savefig(f'{FIG}/fig9_N_P_3d.png', dpi=180); plt.close()

# fig10: lineage distance histogram + host composition
import re
from collections import Counter, defaultdict
alts = defaultdict(dict)
for s_ in spec['proteins']['L']['sites']:
    for alt, lst in s_['alts'].items():
        for a in lst: alts[a['acc']][s_['site']] = alt
sites_all = sorted({s for a in alts.values() for s in a})
dref = {a: sum(1 for k in sites_all if k not in alts[a]) for a in alts}
# careful: alts only stores alt calls; distance = n sites where genome has an alt
dref = {a: len(alts[a]) for a in alts}
fig, axes = plt.subplots(1,2,figsize=(11,3.6))
hmap = {x['accession']: x['host'] for x in spec['genomes']}
byhost = defaultdict(list)
for a,d in dref.items(): byhost[hmap.get(a,'unknown')].append(d)
axes[0].hist([byhost.get('human',[]), byhost.get('bat',[])], bins=25, label=['human','bat'], color=['#C44E52','#55A868'], alpha=0.8)
axes[0].legend(); axes[0].set_xlabel('L variable-site alt count vs reference'); axes[0].set_ylabel('genomes')
axes[0].set_title('Lineage distance from reference by host')
cc = Counter()
for x in spec['genomes']:
    c = str(meta[x['accession']]['country']).split(':')[0]
    cc[(c, x['host'])] += 1
tops = [c for c,_ in Counter({c:sum(v for (c2,_),v in cc.items() if c2==c) for c,_ in cc}).most_common(6)]
widths = {h:[cc.get((c,h),0) for c in tops] for h in ('human','bat','pig','other','unknown')}
bot = [0]*len(tops)
for h,c in zip(('human','bat','pig','other','unknown'),('#C44E52','#55A868','#8172B3','#937860','#999999')):
    axes[1].bar(tops, widths[h], bottom=bot, label=h, color=c)
    bot=[b+w for b,w in zip(bot,widths[h])]
axes[1].legend(fontsize=7); axes[1].set_title('Genomes by country and host'); axes[1].tick_params(axis='x', rotation=30)
plt.tight_layout(); plt.savefig(f'{FIG}/fig10_lineage_geo.png', dpi=200); plt.close()
print('v2 figures done')
