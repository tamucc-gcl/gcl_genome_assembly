"""Current-window placement on native primary graph paths; never infer read support from assembly continuity."""
from collections import defaultdict
import csv
import hashlib
from pathlib import Path
from chimera_controls import peer_relationship


def measure_graph(path, out, queries, intervals, run, threads):
    graph=Path(path);native=out/'native_contigs.fa';segments={};links=[];missing=0
    with graph.open() as handle,native.open('w') as dest:
        for line in handle:
            f=line.rstrip().split('\t')
            if f[0]=='S':
                if len(f)<3 or f[1] in segments:raise ValueError('Invalid/duplicate native GFA segment')
                tags=f[3:];sequence=f[2]
                length=len(sequence) if sequence!='*' else next((int(t[5:]) for t in tags if t.startswith('LN:i:')),0)
                segments[f[1]]=dict(length=length,sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest() if sequence!='*' else None,tags=tags)
                if sequence=='*':missing+=1
                else:dest.write('>'+f[1]+'\n'+sequence+'\n')
            elif f[0]=='L':
                if len(f)<6:raise ValueError('Invalid native GFA link')
                links.append(f)
    if not segments or missing==len(segments):
        native.unlink()
        return {key:dict(graph_status='native_segment_sequences_unavailable',graph_contradiction=None) for key in intervals}
    if sum(r['length'] for r in segments.values())>=8000000000:raise ValueError('Native graph exceeds single-index budget')
    paf=out/'native_paths.paf'
    run(['minimap2','-x','asm5','-I','8G','-c','--eqx','--secondary=yes','-N','50','-p','0.5','-t',threads,native,queries],paf)
    grouped=defaultdict(list)
    for line in paf.read_text().splitlines():
        f=line.split('\t');grouped[f[0]].append(f)
    selected=set();rows=[];result={}
    for key,interval in intervals.items():
        relationship,left,right=peer_relationship(grouped[key],interval['lo']-interval['start'],interval['hi']-interval['start'],anchor=1000,minimum_target=1)
        for hit in (left,right):
            if hit:selected.add(hit[5])
        direct=bool(left and right and any({edge[1],edge[3]}=={left[5],right[5]} for edge in links))
        complete=left is not None and right is not None and not missing
        continuity=relationship=='continuous_context'
        # A native assembled path is a conflicting hypothesis, not independent molecule proof.
        same_segment=bool(left and right and left[5]==right[5])
        contradictory=True if continuity and interval.get('hifi',{}).get('spanning',0)>=2 else None if same_segment or direct or not complete else False
        result[key]=dict(graph_status='screened_primary_contig_paths' if complete else 'native_anchor_placement_unresolved',
            graph_contradiction=contradictory,native_continuity=continuity,direct_native_link=direct,
            native_left=left[5] if left else '.',native_right=right[5] if right else '.',
            native_graph_sha256=None,interpretation='Primary path continuity alone is not read support or biological fusion confirmation; unmeasured unitig paths remain a limitation.')
        rows.append(dict(candidate=key,relationship=relationship,left_segment=result[key]['native_left'],right_segment=result[key]['native_right'],
                         direct_link=direct,graph_status=result[key]['graph_status'],graph_contradiction=contradictory))
    neighbors=set(selected)
    for edge in links:
        if edge[1] in selected or edge[3] in selected:neighbors.update((edge[1],edge[3]))
    with (out/'native_neighborhood.gfa').open('w') as dest:
        dest.write('H\tVN:Z:1.0\n')
        for name in sorted(neighbors):
            if name in segments:
                info=segments[name];tags=[t for t in info['tags'] if not t.startswith('LN:i:')]
                dest.write('\t'.join(['S',name,'*','LN:i:'+str(info['length'])]+tags)+'\n')
        with graph.open() as handle:
            for line in handle:
                f=line.rstrip().split('\t')
                if f[0]=='L' and f[1] in neighbors and f[3] in neighbors:dest.write(line)
                elif f[0]=='A' and len(f)>1 and f[1] in selected:dest.write(line)
    with (out/'native_path_assays.tsv').open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    with graph.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
    for value in result.values():value['native_graph_sha256']=digest
    native.unlink() # Native input remains retained upstream; keep PAF and sequence-free neighborhood.
    return result
