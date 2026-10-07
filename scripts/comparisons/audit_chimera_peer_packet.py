"""Reassess archived peer PAFs without remapping or granting cut authorization."""
import argparse,csv,io,json,sys,tarfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'py_scripts'))
from chimera_controls import peer_anchor_trials

def main():
    p=argparse.ArgumentParser();p.add_argument('--packet',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True);rows=[];summary=[]
    with tarfile.open(a.packet) as archive:
        files={m.name:m for m in archive.getmembers() if m.isfile()}
        for name in files:
            if not name.endswith('control_registry.json'):continue
            root=name.rsplit('/',1)[0]+'/';registry=json.load(archive.extractfile(files[name]))
            provenance=json.load(archive.extractfile(files[root+'provenance.json']))
            assembly=root.rstrip('/').split('/')[-1].replace('.sequence_context','')
            for index,peer in enumerate(provenance['peers']):
                path=root+'peer_%d.paf'%index
                if path not in files:raise ValueError('Missing archived PAF '+path)
                hits={}
                for line in archive.extractfile(files[path]).read().decode().splitlines():
                    f=line.split('\t');hits.setdefault(f[0],[]).append(f)
                for key,interval in registry['intervals'].items():
                    selected,trials=peer_anchor_trials(hits.get(key,[]),interval['lo']-interval['start'],interval['hi']-interval['start'],interval['end']-interval['start'])
                    rows.extend(dict(assembly=assembly,candidate=key,peer=peer['id'],**trial) for trial in trials)
                    relation,left,right=selected
                    labels=peer.get('chromosome_labels',{})
                    lc=labels.get(left[5]) if left else None;rc=labels.get(right[5]) if right else None
                    if relation=='separate_scaffolds':
                        relation='different_chromosomes' if lc and rc and lc!=rc else 'fragmented_same_chromosome' if lc and lc==rc else 'uninformative_chromosome_identity'
                    if relation=='continuous_context' and not lc:relation='uninformative_chromosome_identity'
                    summary.append(dict(assembly=assembly,candidate=key,peer=peer['id'],sample=peer['sample'],auto_evidence=peer.get('auto_evidence',False),relationship=relation,left_chrom=lc or '.',right_chrom=rc or '.'))
    for filename,values in [('anchor_trials.tsv',rows),('peer_summary.tsv',summary)]:
        if values:
            with (out/filename).open('w',newline='') as handle:
                w=csv.DictWriter(handle,fieldnames=list(values[0]),delimiter='\t');w.writeheader();w.writerows(values)
    print('Audited %d peer/interval comparisons; %d anchor trials'%(len(summary),len(rows)))
if __name__=='__main__':main()
