"""Equivalence tests for indexed chromosome projection and process workers."""
import random
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_tracks import LabelIndex,blocks,chromosome_at,position,tracks
from misassembly_core import discover,read_paf,cohort_labels
from misassembly_discover import project_cohort


def naive_tracks(hits,labels):
    """Pre-optimization chromosome-mode reference, including competing placements."""
    segments=[b for f in hits for b in blocks(f)];events=defaultdict(list)
    for i,b in enumerate(segments):events[b['lo']].append((i,1));events[b['hi']].append((i,-1))
    for b in segments:
        value=labels.get(b['target'])
        if isinstance(value,list):
            for label in value:
                for edge in (label['lo'],label['hi']):
                    q=b['lo']+edge-b['t'] if b['strand']=='+' else b['hi']-(edge-b['t'])
                    if b['lo']<q<b['hi']:events[q]
    active=set();result=[];coordinates=sorted(events)
    for index,lo in enumerate(coordinates[:-1]):
        for i,change in events[lo]:
            if change==1:active.add(i)
            else:active.discard(i)
        hi=coordinates[index+1];chrom=None;status='unaligned_or_filtered'
        qualified=[segments[i] for i in active if segments[i]['mapq']>=20 and segments[i]['identity']>=.9 and segments[i]['target'] in labels]
        if qualified:
            best=max(qualified,key=lambda b:(b['mapq'],b['identity']));mid=(lo+hi)//2
            best_chrom=chromosome_at(labels,best['target'],position(best,mid))
            competing=[segments[i] for i in active if segments[i]['identity']>=.95*best['identity'] and (segments[i]['target']!=best['target'] or segments[i]['strand']!=best['strand'] or abs(position(segments[i],mid)-position(best,mid))>1000)]
            if competing and any(chromosome_at(labels,b['target'],position(b,mid))!=best_chrom for b in competing):status='ambiguous'
            else:
                chrom=best_chrom;status=('chromosome_assigned_locus_ambiguous' if competing else 'assigned') if chrom else 'unlabelled_target_block'
        result.append(dict(lo=lo,hi=hi,chrom=chrom,status=status))
    return result


def compact_reference(track):
    out=[]
    for b in track:
        if out and out[-1]['hi']==b['lo'] and (out[-1]['chrom'],out[-1]['status'])==(b['chrom'],b['status']):out[-1]['hi']=b['hi']
        else:out.append(dict(b))
    return out


def hit(target='t',strand='+',cigar='10=1X9=2I10=3D8=',mq=60):
    # 40 query bases, 41 target bases, 95% identity.
    return ['q','100','10','50',strand,target,'100','20','61','38','40',str(mq),'cg:Z:'+cigar]


class ProjectionPerformanceTests(unittest.TestCase):
    def test_index_preserves_conflicts_gaps_and_half_open_edges(self):
        labels={'t':[dict(lo=0,hi=10,chrom='a'),dict(lo=10,hi=20,chrom='a'),dict(lo=15,hi=25,chrom='b'),dict(lo=30,hi=40,chrom='b')]}
        index=LabelIndex(labels)
        for p in range(-1,42):self.assertEqual(index.at('t',p),chromosome_at(labels,'t',p))
        self.assertEqual(index.runs['t'],[(0,15,'a'),(20,25,'b'),(30,40,'b')])

    def test_randomized_tracks_match_previous_algorithm(self):
        rng=random.Random(104)
        for _ in range(80):
            labels={name:[dict(lo=lo,hi=lo+rng.randrange(1,18),chrom=rng.choice(['chr1','chr2'])) for lo in range(0,100,5)] for name in ('t','u')}
            hits=[hit(name,rng.choice(['+','-']),mq=rng.choice([0,20,60])) for name in ('t','u')]
            expected=naive_tracks(hits,labels)
            self.assertEqual(tracks(hits,labels,mode='chromosome',min_mapq=20,minimum_target=0),expected)
            actual=tracks(hits,labels,mode='chromosome',min_mapq=20,minimum_target=0,compact=True)
            self.assertEqual(actual,compact_reference(expected))
            peer=dict(id='peer',sample='other',eligible=True)
            self.assertEqual(discover([(peer,{'q':expected})],'focal',minimum=2),discover([(peer,{'q':actual})],'focal',minimum=2))

    def test_compact_cigar_keeps_indel_boundaries_on_both_strands(self):
        for strand in ('+','-'):
            original=blocks(hit(strand=strand));merged=blocks(hit(strand=strand),compact=True)
            self.assertEqual(len(merged),3)
            for p in range(10,50):
                left=[position(b,p) for b in original if b['lo']<=p<b['hi']]
                right=[position(b,p) for b in merged if b['lo']<=p<b['hi']]
                self.assertEqual(left,right)

    def test_filter_before_reverse_cigar_parsing(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.paf';p.write_text('\t'.join(hit())+'\n'+'\t'.join(hit(target='outside',cigar='malformed'))+'\n',encoding='utf-8')
            self.assertEqual(len(list(read_paf(p,True,{'t'}))),1)

    def test_catalog_merges_only_equivalent_supporting_individuals(self):
        with tempfile.TemporaryDirectory() as d:
            catalog={'focal':dict(id='focal',sample='focal',taxid='1',scope=['q'],scope_lengths={'q':100},eligible=True,chromosome_labels={})}
            pairs=[]
            for i in range(3):
                name='peer'+str(i);p=Path(d)/(name+'.paf');records=[]
                for start in ([10,30] if i<2 else [30]):
                    f=['q','100',str(start),str(start+20),'+','t','100',str(start),str(start+20),'20','20','60','cg:Z:20=']
                    records.append('\t'.join(f))
                p.write_text('\n'.join(records)+'\n',encoding='utf-8')
                catalog[name]=dict(id=name,sample=name,taxid='1',scope=['t'],eligible=True,chromosome_labels={'t':'chr1'})
                pairs.append(dict(a='focal',b=name,path=str(p)))
            result=cohort_labels(catalog,pairs)['focal']['chromosome_labels']['q']
            self.assertEqual([(b['lo'],b['hi'],b['source_individuals']) for b in result],[(10,30,['peer0','peer1']),(30,50,['peer0','peer1','peer2'])])

    def test_parallel_peer_projection_matches_serial(self):
        with tempfile.TemporaryDirectory() as d:
            pairs=[];catalog={'focal':dict(id='focal',sample='focal',taxid='1',scope=['q'],scope_lengths={'q':100},chromosome_labels={})}
            for i in range(3):
                name='peer'+str(i);p=Path(d)/(name+'.paf');p.write_text('\t'.join(hit())+'\n',encoding='utf-8')
                catalog[name]=dict(id=name,sample=name,taxid='1',eligible=True,chromosome_labels={'t':'chr1'})
                pairs.append(dict(a='focal',b=name,path=str(p)))
            self.assertEqual(project_cohort(catalog,'focal',pairs,1),project_cohort(catalog,'focal',pairs,2))


if __name__=='__main__':unittest.main()
