#!/usr/bin/env python3
"""G2: per-site amino-acid mutation spectrum of NiV L/P/N across all genomes
vs reference NC_002728.1, with accession-level provenance and host
stratification. Also checks catalytic-residue invariance (G1 follow-through).
Genomes lacking CDS annotation: L/P/N recovered by 6-frame translation +
local alignment to reference proteins (provenance flagged 'orf_recovered').
"""
import json, os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RES = os.path.join(HERE, '..', 'results')
REF = 'NC_002728.1'

def classify(f):
    q = f.qualifiers
    gene = (q.get('gene') or [''])[0].upper()
    prod = (q.get('product') or [''])[0].lower()
    if gene == 'N' or 'nucleocapsid' in prod or 'nucleoprotein' in prod: return 'N'
    if gene in ('P','P/V','P/V/C') or 'phosp' in prod or prod == 'p protein': return 'P'
    if gene == 'L' or 'polymerase' in prod or prod == 'l protein': return 'L'
    return None

def host_group(h):
    h = h.lower()
    if 'homo sapiens' in h: return 'human'
    if 'pteropus' in h or h == 'bat': return 'bat'
    if 'sus scrofa' in h: return 'pig'
    if h == 'unknown': return 'unknown'
    return 'other'

aln = PairwiseAligner()
aln.mode = 'global'
aln.match_score = 2; aln.mismatch_score = -1
aln.open_gap_score = -5; aln.extend_gap_score = -1

def spectrum_pair(ref, qry):
    """Return dict site(1-based ref) -> alt aa for mismatches, plus identity."""
    a = aln.align(ref, qry)[0]
    coords = a.aligned
    subs = {}
    rpos = 0; qpos = 0
    rseq = a[0]; qseq = a[1]
    matches = 0; aligned_len = 0
    for i in range(len(rseq)):
        r, q = rseq[i], qseq[i]
        if r != '-':
            rpos += 1
        if q != '-':
            qpos += 1
        if r != '-' and q != '-':
            aligned_len += 1
            if r == q:
                matches += 1
            elif q in 'XBZ':
                subs[rpos] = q  # tracked separately, not a real substitution
            else:
                subs[rpos] = q
    ident = matches / aligned_len if aligned_len else 0.0
    return subs, ident

def recover_orf(genome, refprot):
    """6-frame scan: find best local alignment of refprot in genome."""
    la = PairwiseAligner(); la.mode = 'local'
    la.match_score = 2; la.mismatch_score = -1
    la.open_gap_score = -10; la.extend_gap_score = -1
    best = None
    for strand, nuc in ((+1, genome), (-1, genome.reverse_complement())):
        s = str(nuc)
        for frame in range(3):
            prot = str(Seq(s[frame:frame + 3 * ((len(s) - frame) // 3)]).translate())
            a = la.align(prot, refprot)[0]
            if best is None or a.score > best[0]:
                seg = a.aligned[0]
                best = (a.score, prot[seg[0][0]:seg[0][1]], a)
    score, qseq, a = best
    ident = sum(1 for x, y in zip(a[0], a[1]) if x == y) / max(1, len(a[0]))
    return qseq, ident

def main():
    records = {r.id: r for r in SeqIO.parse(os.path.join(DATA, 'genomes', 'niv_genomes.gb'), 'genbank')}
    manifest = json.load(open(os.path.join(DATA, 'manifest.json')))
    hosts = {g['accession']: host_group(g['host']) for g in manifest['genomes']}
    ref = records[REF]
    refprot = {}
    for f in ref.features:
        if f.type == 'CDS':
            g = classify(f)
            tr = (f.qualifiers.get('translation') or [''])[0]
            if g and tr and g not in refprot:
                refprot[g] = tr
    out = {'reference': REF, 'proteins': {}, 'genomes': [], 'excluded': []}
    spec = {g: {} for g in 'LPN'}  # spec[g][site] = {'ref':aa,'alts':{alt:[{acc,host}]}}
    for acc, rec in sorted(records.items()):
        if acc == REF: continue
        hg = hosts.get(acc, 'unknown')
        prots = {}
        for f in rec.features:
            if f.type == 'CDS':
                g = classify(f)
                tr = (f.qualifiers.get('translation') or [''])[0]
                if g and tr and g not in prots:
                    prots[g] = (tr, 'cds')
        missing = [g for g in 'LPN' if g not in prots]
        prov = 'cds'
        if missing:
            # ORF recovery fallback for any protein missing from annotation
            for g in missing[:]:
                q, ident = recover_orf(rec.seq, refprot[g])
                if ident >= 0.90 and len(q) >= 0.95 * len(refprot[g]):
                    prots[g] = (q, 'orf_recovered')
                    missing.remove(g)
            if any(v[1] == 'orf_recovered' for v in prots.values()):
                prov = 'mixed_orf'
        if missing:
            out['excluded'].append({'accession': acc, 'host': hg, 'missing': missing,
                                    'n_cds': len([f for f in rec.features if f.type == 'CDS'])})
            continue
        ginfo = {'accession': acc, 'host': hg, 'provenance': prov, 'identity': {}}
        for g in 'LPN':
            qseq, how = prots[g]
            subs, ident = spectrum_pair(refprot[g], qseq)
            ambig = {k: v for k, v in subs.items() if v in 'XBZ'}
            subs = {k: v for k, v in subs.items() if v not in 'XBZ'}
            ginfo['identity'][g] = round(ident, 5)
            ginfo['n_subs_' + g] = len(subs)
            if ambig:
                ginfo['n_ambiguous_' + g] = len(ambig)
            for site, alt in subs.items():
                s = spec[g].setdefault(site, {'ref': refprot[g][site-1], 'alts': {}})
                s['alts'].setdefault(alt, []).append({'acc': acc, 'host': hg})
        out['genomes'].append(ginfo)
        print(acc, hg, prov, {g: ginfo['n_subs_' + g] for g in 'LPN'})
    # finalize spectrum
    for g in 'LPN':
        sites = []
        for site, d in spec[g].items():
            n = sum(len(v) for v in d['alts'].values())
            hosts_alt = sorted({a['host'] for v in d['alts'].values() for a in v})
            sites.append({'site': site, 'ref': d['ref'], 'n_carriers': n,
                          'host_groups': hosts_alt, 'alts': d['alts']})
        sites.sort(key=lambda x: x['site'])
        out['proteins'][g] = {'length': len(refprot[g]), 'n_variable_sites': len(sites), 'sites': sites}
    # catalytic invariance (G1 follow-through)
    cat = {'L': [722, 831, 832, 833, 834, 1165, 1821, 1940, 1976, 2013]}
    inv = {}
    for g, positions in cat.items():
        inv[g] = {}
        for pos in positions:
            d = spec[g].get(pos)
            inv[g][str(pos)] = ('INVARIANT' if d is None else
                                'VARIABLE: ' + json.dumps(d['alts']))
    out['catalytic_invariance'] = inv
    # host-exclusive sites
    strat = {}
    for g in 'LPN':
        rows = []
        for s in out['proteins'][g]['sites']:
            hgset = set(s['host_groups'])
            known = hgset - {'unknown', 'other'}
            if known and len(known) == 1 and s['n_carriers'] >= 2:
                rows.append({'site': s['site'], 'ref': s['ref'],
                             'exclusive_host': list(known)[0],
                             'n_carriers': s['n_carriers'],
                             'alts': {k: len(v) for k, v in s['alts'].items()}})
        strat[g] = rows
    out['host_exclusive_sites'] = strat
    json.dump(out, open(os.path.join(RES, 'g2_spectrum.json'), 'w'), indent=1)
    print('=== G2 summary ===')
    print('analyzed genomes:', len(out['genomes']), 'excluded:', len(out['excluded']))
    for g in 'LPN':
        print(g, 'len', out['proteins'][g]['length'],
              'variable sites:', out['proteins'][g]['n_variable_sites'],
              'host-exclusive sites (>=2 carriers):', len(strat[g]))
    print('catalytic invariance:', json.dumps(inv))

if __name__ == '__main__':
    main()
