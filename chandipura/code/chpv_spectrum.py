#!/usr/bin/env python3
"""G2: per-site mutation spectrum with accession-level provenance + host stratification."""
import json, os, math
from collections import Counter, defaultdict
from scipy.stats import fisher_exact
from chpv_profile import (ROOT, REF, HUMAN, VECTOR, HEDGEHOG, LAB, EXCLUDED_DUP,
                          load_proteins, group_of, build_alignment, entropy, DOMAINS_L, LANDMARKS)

prot = load_proteins()
ALL = [a for a in prot['L'] if a not in EXCLUDED_DUP]
NATURAL = [a for a in ALL if group_of(a) in ('human', 'vector', 'hedgehog')]
print(f'genomes in spectrum: {len(ALL)} (27 = 28 - 1 RefSeq dup); natural isolates: {len(NATURAL)}')
print('human:', len(HUMAN), 'vector:', len(VECTOR), 'hedgehog:', len(HEDGEHOG), 'lab:', len(LAB), '+reference')

aln = {g: build_alignment(g, prot, ALL) for g in ['L', 'P', 'N']}
refseq = {g: prot[g][REF] for g in ['L', 'P', 'N']}

spectrum = {}   # gene -> list of site dicts
subs_table = [] # every substitution with provenance
for gene in ['N', 'P', 'L']:
    ref = refseq[gene]
    L = len(ref)
    sites = []
    for i in range(L):
        col = {acc: aln[gene][acc][i] for acc in ALL}
        counts = Counter(col.values())
        h = entropy(counts)
        var_accs = {acc: aa for acc, aa in col.items() if aa != ref[i]}
        site = dict(pos=i+1, ref=ref[i], entropy=round(h, 4),
                    n_variant=len(var_accs), groups=Counter(group_of(a) for a in var_accs))
        sites.append(site)
        for acc, aa in var_accs.items():
            subs_table.append(dict(gene=gene, pos=i+1, ref=ref[i], alt=aa,
                                   accession=acc, group=group_of(acc)))
    spectrum[gene] = sites
    n_var_sites = sum(1 for s in sites if s['n_variant'] > 0)
    print(f'{gene}: {L} sites, {n_var_sites} variable ({n_var_sites/L*100:.1f}%), '
          f'mean entropy {sum(s["entropy"] for s in sites)/L:.4f}')

# verified counts
tot_subs = len(subs_table)
per_gene = Counter(s['gene'] for s in subs_table)
print('total substitutions (all genomes vs I653514):', tot_subs, dict(per_gene))
check = sum(s['n_variant'] for g in spectrum for s in spectrum[g])
print('reconciliation: table rows', tot_subs, '== site sums', check, '->', tot_subs == check)

# ---- host-stratified: human vs vector, natural isolates only, per-site Fisher ----
host_hits = []
for gene in ['N', 'P', 'L']:
    ref = refseq[gene]
    for i in range(len(ref)):
        hu = sum(1 for a in HUMAN if aln[gene][a][i] != ref[i])
        ve = sum(1 for a in VECTOR if aln[gene][a][i] != ref[i])
        if hu == 0 and ve == 0: continue
        _, p = fisher_exact([[hu, len(HUMAN)-hu], [ve, len(VECTOR)-ve]])
        if hu > 0 or ve > 0:
            host_hits.append(dict(gene=gene, pos=i+1, ref=ref[i],
                                  human_var=hu, vector_var=ve, fisher_p=round(p, 4)))
sig = [h for h in host_hits if h['fisher_p'] < 0.05]
print(f'\nvariable sites in human/vector comparison: {len(host_hits)}; Fisher p<0.05: {len(sig)}')
for h in sorted(sig, key=lambda x: x['fisher_p'])[:25]:
    print(f"  {h['gene']}{h['pos']} {h['ref']}: human {h['human_var']}/5 vector {h['vector_var']}/17 p={h['fisher_p']}")

# human-specific (variant in ALL 5 human, never in vector) and vector-specific
def group_specific(gene, group_accs, other_accs):
    ref = refseq[gene]
    out = []
    for i in range(len(ref)):
        if all(aln[gene][a][i] != ref[i] and aln[gene][a][i] == aln[gene][group_accs[0]][i] for a in group_accs):
            if all(aln[gene][a][i] == ref[i] for a in other_accs):
                out.append((i+1, ref[i], aln[gene][group_accs[0]][i]))
    return out
hu_spec = {g: group_specific(g, HUMAN, VECTOR) for g in ['N','P','L']}
ve_spec = {g: group_specific(g, VECTOR, HUMAN) for g in ['N','P','L']}
for g in ['N','P','L']:
    print(f'{g}: human-uniform-specific sites {len(hu_spec[g])}, vector-uniform-specific sites {len(ve_spec[g])}')

# ---- domain-level stats for L ----
print('\nL domain-level mean entropy:')
dom_stats = []
assigned = set()
for name, lo, hi in DOMAINS_L:
    vals = [spectrum['L'][i-1]['entropy'] for i in range(lo, hi+1)]
    dom_stats.append(dict(domain=name, start=lo, end=hi,
                          mean_entropy=round(sum(vals)/len(vals), 4),
                          variable_sites=sum(1 for v in vals if v > 0)))
    assigned.update(range(lo, hi+1))
    print(f'  {name:10s} {lo:4d}-{hi:4d}: meanH={dom_stats[-1]["mean_entropy"]:.4f} var={dom_stats[-1]["variable_sites"]}/{hi-lo+1}')
rest = [spectrum['L'][i-1]['entropy'] for i in range(1, len(refseq['L'])+1) if i not in assigned]
print(f'  outside domains: meanH={sum(rest)/len(rest):.4f} var={sum(1 for v in rest if v>0)}/{len(rest)}')

json.dump(dict(spectrum=spectrum, subs_table=subs_table, host_hits=host_hits,
               human_specific={g: hu_spec[g] for g in hu_spec},
               vector_specific={g: ve_spec[g] for g in ve_spec},
               l_domain_stats=dom_stats,
               counts=dict(total_subs=tot_subs, per_gene=dict(per_gene),
                           reconciled=(tot_subs == check))),
         open(os.path.join(ROOT, 'results/g2_spectrum.json'), 'w'))
print('\nwrote results/g2_spectrum.json')
