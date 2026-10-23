"""Production chromosome-boundary ownership and query-lane filtering utilities.

Promoted unchanged from the validated diagnostic implementation. No assembly or
chromosome names are built into this policy.
"""
from array import array
from collections import Counter
import gzip
import hashlib
from pangenome_source_audit import read_table, write_table


def lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    return opener(path,'rt',encoding='utf-8')


def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def node_owners(gfa, output):
    """Same chromosome membership as reference-boundary BFS, without per-node BFS.

    Rank-zero nodes are boundaries. Only non-reference nodes are unioned; reference
    nodes must NOT join adjacent chromosomes through a composite graph edge.
    Components bordering multiple chromosomes are explicit errors when used.
    """
    ids, names, chromosomes = {}, [], {}
    parent, rank, ref = array('i'), array('B'), array('i')
    with lines(gfa) as handle:
        for line in handle:
            if not line.startswith('S\t'): continue
            f = line.rstrip().split('\t'); name = f[1]
            if name in ids: raise ValueError('Duplicate graph node: '+name)
            idx = len(parent); ids[name] = idx; names.append(name)
            tags = {t[:2]:t[5:] for t in f[3:] if len(t)>5 and t[2]==':' and t[4]==':'}
            chrom = tags.get('SN') if tags.get('SR')=='0' else None
            if chrom: chromosomes.setdefault(chrom,len(chromosomes)+1)
            parent.append(idx); rank.append(0); ref.append(chromosomes[chrom] if chrom else 0)
    if not chromosomes: raise ValueError('SV rGFA lacks rank-zero reference chromosome tags')
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    def edges():
        with lines(gfa) as handle:
            for line in handle:
                if line.startswith('L\t'):
                    f = line.split('\t')
                    if f[1] not in ids or f[3] not in ids: raise ValueError('Unknown graph link node')
                    yield ids[f[1]],ids[f[3]]
    for a,b in edges():
        if ref[a] or ref[b]: continue
        a,b = root(a),root(b)
        if a==b: continue
        if rank[a]<rank[b]: a,b=b,a
        parent[b]=a
        if rank[a]==rank[b]: rank[a]+=1
    owner = array('i',[0])*len(parent)
    for a,b in edges():
        if bool(ref[a]) == bool(ref[b]): continue
        if ref[b]: a,b=b,a
        component = root(b); value = owner[component]
        owner[component] = ref[a] if value in (0,ref[a]) else -1
    labels = {value:name for name,value in chromosomes.items()}
    def rows():
        for i,name in enumerate(names):
            c = ref[i] or owner[root(i)]
            yield dict(node=name,reference_chromosome=labels.get(c,''),
                       status='anchored' if c>0 else 'conflicting' if c<0 else 'unanchored')
    write_table(output,rows(),['node','reference_chromosome','status'])
    return dict(nodes=len(parent),reference_chromosomes=len(chromosomes))


def lane_paf(source,owners_path,dest,ledger):
    owners={r['node']:(r['reference_chromosome'],r['status']) for r in read_table(owners_path)}
    lanes={};counts=Counter();matches=Counter()
    with lines(source) as src,open(dest,'w',encoding='utf-8') as out:
        for line in src:
            f=line.rstrip('\n').split('\t')
            if len(f)<12: raise ValueError('Malformed lane PAF')
            node=f[5].split('|',1)[1] if f[5].startswith('id=_MINIGRAPH_|') else f[5]
            if node not in owners or owners[node][1]!='anchored':
                raise ValueError('Cannot safely assign mapping node to one chromosome: '+node)
            key=f[0],owners[node][0]
            if key not in lanes: lanes[key]='gcl_lane_'+str(len(lanes)+1)
            f[0]=lanes[key];counts[key]+=1;matches[key]+=int(f[9])
            out.write('\t'.join(f)+'\n')
    write_table(ledger,(dict(lane=lane,query=key[0],reference_chromosome=key[1],
        records=counts[key],matching_bases=matches[key]) for key,lane in lanes.items()),
        ['lane','query','reference_chromosome','records','matching_bases'])
    return {lane:key[0] for key,lane in lanes.items()}


def restore_names(source,dest,names):
    with open(source,encoding='utf-8') as src,open(dest,'w',encoding='utf-8') as dst:
        for line in src:
            query,rest=line.split('\t',1)
            if query not in names: raise ValueError('Unknown filtering lane: '+query)
            dst.write(names[query]+'\t'+rest)
