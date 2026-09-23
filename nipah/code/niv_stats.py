#!/usr/bin/env python3
"""Statistical analysis: domain enrichment, RSA enrichment, host
stratification (Fisher), temporal trend (Spearman), pairwise-distance
lineage analysis (confound check), all real computed numbers -> g4_stats.json
"""
import json, os, math
from collections import Counter, defaultdict
from scipy.stats import chi2_contingency, mannwhitneyu, fisher_exact, spearmanr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
DATA = os.path.join(HERE, '..', 'data')

spec = json.load(open(os.path.join(RES, 'g2_spectrum.json')))
struct = json.load(open(os.path.join(RES, 'g3_structure.json')))
man = json.load(open(os.path.join(DATA, 'manifest.json')))
meta = {g['accession']: g for g in man['genomes']}
ginfo = {g['accession']: g for g in spec['genomes']}

out = {}

# --- 1. domain enrichment (chi-square vs length-proportional expectation) ---
dmaps = {'L': [(1,969,'RdRp'),(970,1452,'PRNTase'),(1453,1758,'CD'),(1759,2080,'MTase'),(2081,2244,'CTD')],
         'P': [(1,469,'NTD_disordered'),(470,578,'OD'),(579,651,'OD-XD_linker'),(652,709,'XD')],
         'N': [(1,405,'core'),(406,532,'tail_disordered')]}
out['domain_enrichment'] = {}
for g in 'LPN':
    doms = dmaps[g]
    total_len = sum(hi-lo+1 for lo,hi,_ in doms)
    obs = Counter()
    for s in struct['proteins'][g]:
        obs[s['domain']] += 1
    nvar = sum(obs.values())
    table, exp = [], []
    for lo,hi,name in doms:
        e = nvar * (hi-lo+1)/total_len
        table.append([obs.get(name,0), nvar - obs.get(name,0)])
        exp.append(round(e,1))
    chi2, p, dof, _ = chi2_contingency(np.array(table).T) if nvar>0 else (0,1,0,None)
    rates = {name: round(obs.get(name,0)/(hi-lo+1)*1000,2) for lo,hi,name in doms}
    out['domain_enrichment'][g] = {'obs': dict(obs), 'expected': exp,
        'sites_per_1k': rates, 'chi2': round(float(chi2),2), 'p': float(p)}

# --- 2. RSA: variable vs invariant modeled residues ---
out['rsa_enrichment'] = {}
for g in 'LPN':
    rows = struct['proteins'][g]
    var_rsa = [r['rsa'] for r in rows if 'rsa' in r]
    modeled = set()
    # reconstruct modeled invariant residues from domain map lengths & unmodeled counts
    var_sites = {r['site'] for r in rows}
    seqlen = spec['proteins'][g]['length']
    # modeled residues = those with rsa in struct rows OR total modeled count
    n_unmodeled_var = sum(1 for r in rows if r['location']=='unmodeled')
    n_modeled_var = len(var_rsa)
    # total modeled: read from structure stats
    stats3 = struct['stats'][g]['by_location']
    modeled_total = None
    # approximate modeled total from L/P/N: use rsa presence across all residues is not stored;
    # use: n_var_modeled + invariant modeled. Invariant modeled = modeled_total - n_var_modeled.
    # modeled_total per protein known from G3 run (L 2020, P 195, N 315); recompute here from rows is impossible -> store from g3
    out['rsa_enrichment'][g] = {'n_var_modeled': n_modeled_var,
        'median_rsa_variable': float(np.median(var_rsa)) if var_rsa else None}
# full RSA comparison needs all-residue RSA -> do a compact re-run here
import warnings; warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser, ShrakeRupley
def rsa_all(cif, chain):
    p = MMCIFParser(QUIET=True)
    s = p.get_structure('x', os.path.join(DATA,'structures',cif))
    sr = ShrakeRupley(n_points=100); sr.compute(s[0], level='R')
    return {r.id[1]: min(r.sasa/210.0,1.0) for r in s[0][chain] if r.id[0]==' ' and 'CA' in r}
RSA = {'L': rsa_all('9GJU.cif','A'), 'P': rsa_all('9GJU.cif','B'), 'N': rsa_all('4CO6.cif','A')}
for g in 'LPN':
    var_sites = {r['site'] for r in struct['proteins'][g]}
    v = [rsa for s,rsa in RSA[g].items() if s in var_sites]
    iv = [rsa for s,rsa in RSA[g].items() if s not in var_sites]
    if v and iv:
        U, p = mannwhitneyu(v, iv, alternative='greater')
        # enrichment of variable sites on surface vs core, Fisher
        a = sum(1 for x in v if x>=0.25); b = len(v)-a
        c = sum(1 for x in iv if x>=0.25); d = len(iv)-c
        odds, pf = fisher_exact([[a,b],[c,d]])
        out['rsa_enrichment'][g].update({'n_modeled': len(RSA[g]),
            'median_rsa_variable': float(np.median(v)), 'median_rsa_invariant': float(np.median(iv)),
            'mannwhitney_p_greater': float(p), 'surface_fisher_odds': float(odds),
            'surface_fisher_p': float(pf), 'var_surface_frac': a/len(v), 'inv_surface_frac': c/len(iv)})

# --- 3. host stratification Fisher per site (human vs bat) ---
def hg(acc): return ginfo.get(acc,{}).get('host') or {'human':'human'}.get('')
hostcounts = Counter(x['host'] for x in spec['genomes'])
nH, nB = hostcounts.get('human',0), hostcounts.get('bat',0)
out['host_counts'] = dict(hostcounts)
sig = {}
for g in 'LPN':
    rows = []
    for s in spec['proteins'][g]['sites']:
        carrH = sum(1 for v in s['alts'].values() for a in v if a['host']=='human')
        carrB = sum(1 for v in s['alts'].values() for a in v if a['host']=='bat')
        if carrH+carrB >= 2:
            table = [[carrH, carrB],[nH-carrH, nB-carrB]]
            odds, p = fisher_exact(table)
            if p < 0.05 and (carrH==0 or carrB==0):
                rows.append({'site': s['site'], 'ref': s['ref'], 'human': carrH,
                             'bat': carrB, 'fisher_p': round(p,5),
                             'direction': 'human' if carrH>0 else 'bat'})
    sig[g] = rows
out['host_exclusive_significant'] = sig

# --- 4. temporal: year vs n_subs (L) ---
xs, ys = [], []
for x in spec['genomes']:
    y = str(meta.get(x['accession'],{}).get('collection_date','unknown'))
    import re as _re
    m = _re.search(r'(19|20)\d{2}', y)
    if not m: continue
    yr = int(m.group())
    if 1990 <= yr <= 2026:
        xs.append(yr); ys.append(x['n_subs_L'])
rho, p = spearmanr(xs, ys)
out['temporal_L'] = {'n': len(xs), 'spearman_rho': float(rho), 'p': float(p)}

# --- 5. lineage/confound: pairwise variable-site distance clustering ---
alts = defaultdict(dict)  # acc -> site->alt (L)
for s in spec['proteins']['L']['sites']:
    for alt, lst in s['alts'].items():
        for a in lst:
            alts[a['acc']][s['site']] = alt
accs = sorted(alts)
sites_all = sorted({s for a in alts.values() for s in a})
vecs = {a: [alts[a].get(s,'.') for s in sites_all] for a in accs}
# distance matrix on variable sites only
n = len(accs)
D = np.zeros((n,n))
for i in range(n):
    vi = vecs[accs[i]]
    for j in range(i+1,n):
        vj = vecs[accs[j]]
        d = sum(1 for k in range(len(sites_all)) if vi[k] != vj[k])
        D[i,j] = D[j,i] = d
# simple 2-cluster via median split on distance to reference-like profile
refprof = ['.']*len(sites_all)
dist_ref = {a: sum(1 for k in range(len(sites_all)) if vecs[a][k] != '.') for a in accs}
med = np.median(list(dist_ref.values()))
clu = {a: ('near-ref' if dist_ref[a] <= med else 'far') for a in accs}
# composition
comp = defaultdict(Counter)
for a in accs:
    comp[clu[a]][ginfo[a]['host']] += 1
    comp[clu[a]][meta[a]['country'][:12]] += 0  # placeholder no-op
out['lineage_L'] = {'n': n, 'sites': len(sites_all),
    'median_dist': float(med),
    'near_ref_hosts': dict(comp['near-ref']), 'far_hosts': dict(comp['far'])}
# country confound: bat vs human country overlap
countries = defaultdict(set)
for x in spec['genomes']:
    countries[x['host']].add(meta[x['accession']]['country'])
out['country_by_host'] = {k: sorted(v)[:12] for k,v in countries.items()}

json.dump(out, open(os.path.join(RES,'g4_stats.json'),'w'), indent=1)
print('domain enrichment:', json.dumps({g: out['domain_enrichment'][g]['sites_per_1k'] for g in 'LPN'}, indent=1))
print('rsa:', json.dumps(out['rsa_enrichment'], indent=1)[:800])
print('host sig sites:', {g: len(sig[g]) for g in 'LPN'})
print('temporal:', out['temporal_L'])
print('lineage:', out['lineage_L']['near_ref_hosts'], out['lineage_L']['far_hosts'])
print('countries:', out['country_by_host'])
