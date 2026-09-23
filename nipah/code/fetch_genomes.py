#!/usr/bin/env python3
"""G0 byte-lock: fetch all complete Nipah virus genomes from NCBI eutils and
build per-accession FASTAs + md5 manifest. Re-run to regenerate data/."""
import urllib.request, urllib.parse, time, json, hashlib, io, os
from Bio import SeqIO

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', 'genomes')
PER = os.path.join(DATA, 'per_accession')
TERM = 'txid121791[Organism] AND 17000:19000[SLEN]'
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def get(url, data=None):
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, data=data)
            return urllib.request.urlopen(req, timeout=180).read()
        except Exception as e:
            wait = min(2 ** attempt * 2, 30)
            print('retry', attempt, type(e).__name__, 'sleep', wait)
            time.sleep(wait)
    raise

def main():
    os.makedirs(PER, exist_ok=True)
    url = EUTILS + 'esearch.fcgi?' + urllib.parse.urlencode(
        {'db':'nuccore','term':TERM,'retmax':500,'retmode':'json'})
    ids = json.loads(get(url))['esearchresult']['idlist']
    print('esearch ids:', len(ids))
    gb_path = os.path.join(DATA,'niv_genomes.gb')
    fa_path = os.path.join(DATA,'niv_genomes_raw.fasta')
    if os.path.exists(gb_path) and os.path.exists(fa_path):
        print('using cached bulk files')
        gb_parts = [open(gb_path,'rb').read()]; fa_parts = [open(fa_path,'rb').read()]
    else:
        gb_parts, fa_parts = [], []
        for i in range(0, len(ids), 50):
            chunk = ids[i:i+50]
            for rettype, sink in (('gbwithparts', gb_parts), ('fasta', fa_parts)):
                data = urllib.parse.urlencode({'db':'nuccore','id':','.join(chunk),
                                           'rettype':rettype,'retmode':'text'}).encode()
                sink.append(get(EUTILS+'efetch.fcgi', data))
                time.sleep(0.4)
    open(os.path.join(DATA,'niv_genomes.gb'),'wb').write(b''.join(gb_parts))
    open(os.path.join(DATA,'niv_genomes_raw.fasta'),'wb').write(b''.join(fa_parts))
    manifest = []
    for rec in SeqIO.parse(io.StringIO(b''.join(gb_parts).decode()), 'genbank'):
        seq = str(rec.seq).upper()
        fasta = f'>{rec.id} {rec.description}\n{seq}\n'
        path = os.path.join(PER, rec.id + '.fasta')
        open(path,'w').write(fasta)
        md5 = hashlib.md5(fasta.encode()).hexdigest()
        srcf = next((f for f in rec.features if f.type == 'source'), None)
        src = srcf.qualifiers if srcf else {}
        cds = [f for f in rec.features if f.type=='CDS']
        manifest.append({
            'accession': rec.id,
            'length': len(seq),
            'host': (src.get('host') or ['unknown'])[0],
            'country': (src.get('country') or ['unknown'])[0],
            'collection_date': (src.get('collection_date') or ['unknown'])[0],
            'isolate': (src.get('isolate') or src.get('strain') or ['unknown'])[0],
            'n_cds': len(cds),
            'md5': md5,
            'file': 'data/genomes/per_accession/%s.fasta' % rec.id,
        })
    manifest.sort(key=lambda r: r['accession'])
    out = {'project':'NIV-STRUCT','taxid':'121791','n_genomes':len(manifest),
           'fetched':'2026-09-23','source':'NCBI nuccore eutils efetch',
           'query': TERM, 'genomes': manifest}
    with open(os.path.join(DATA,'..','manifest.json'),'w') as fh:
        json.dump(out, fh, indent=2)
    hosts = {}
    for m in manifest: hosts[m['host']] = hosts.get(m['host'],0)+1
    print('genomes:', len(manifest), 'hosts:', hosts)

if __name__ == '__main__':
    main()
