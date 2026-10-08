import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_tracks import tracks,augment_peer_labels
from chimera_markdown import decision_context
from chimera_review import generate
from chimera_controls import scan_contacts
from chimera_report import families


def paf(query='q',target='composite',strand='+',lo=0,hi=200000):
    return [query,'200000','0','200000',strand,target,'10000000',str(lo),str(hi),'200000','200000','60','cg:Z:200000M']


class ReportRevision(unittest.TestCase):
    def test_composite_target_is_split_into_chromosome_blocks(self):
        labels={'composite':[dict(lo=0,hi=100000,chrom='chr7'),dict(lo=100000,hi=200000,chrom='chr12')]}
        forward=tracks([paf()],labels)
        reverse=tracks([paf(strand='-')],labels)
        self.assertEqual([(x['lo'],x['hi'],x['chrom']) for x in forward],[(0,100000,'chr7'),(100000,200000,'chr12')])
        self.assertEqual([x['chrom'] for x in reverse],['chr12','chr7'])

    def test_ambiguous_reference_block_does_not_become_whole_composite_label(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'ref.paf'
            a=paf(query='composite',target='r7');a[2:4]=['0','200000']
            b=paf(query='composite',target='r12');b[2:4]=['100000','300000'];b[1]='300000'
            path.write_text('\t'.join(a)+'\n'+'\t'.join(b)+'\n',encoding='utf-8')
            peer=dict(reference_paf=str(path),chromosome_labels={},reference_labels={'r7':'chr7','r12':'chr12'})
            augment_peer_labels(peer)
            self.assertEqual(peer['chromosome_labels']['composite'],[dict(lo=0,hi=100000,chrom='chr7'),dict(lo=200000,hi=300000,chrom='chr12')])

    def test_consistent_partial_assignments_are_evidence_not_no_evidence(self):
        side=lambda c:dict(chrom=c,aligned_bp=80000,qualified=False)
        tracks_=[dict(peer=s+'_hap1',sample=s,auto_evidence=True,relationship='uninformative',left=side('chr7'),right=side('chr12')) for s in ('A','B','C')]
        result=decision_context(dict(chromosome_tracks=tracks_),'self')
        self.assertIn('3 independent individuals',result['evidence_for_cut'])
        self.assertIn('0 independent individuals separate',result['local_confirmation'])

    def test_completed_empty_contact_scan_emits_measured_zeros(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'pairs';path.write_text('# header\n')
            counts,totals,audit=scan_contacts(path,{},dict(c=dict(scaffold='s',lo=10,hi=20)),[dict(library_id='A',qname_prefix='A_')])
            self.assertEqual(dict(counts['A','c']),dict(cross=0,left_within=0,right_within=0,left_ends=0,right_ends=0))

    def test_report_is_compact_and_restores_scaffold_figures(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);ctx=root/'context';ctx.mkdir();supp=root/'supp';supp.mkdir()
            record=dict(scaffold='s',assembly_sha256='a'*64,transition_lo='100',transition_hi='200',left_chrom='chr7',right_chrom='chr12')
            calls=root/'calls.tsv'
            with calls.open('w',newline='',encoding='utf-8') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(record),delimiter='\t');writer.writeheader();writer.writerow(record)
            key=hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:20]
            (ctx/'provenance.json').write_text(json.dumps(dict(sample='self',sha256='a'*64,intervals={})))
            (ctx/'decision_measurements.json').write_text(json.dumps({key:dict(native_continuity=True,native_left='n1',native_right='n1',large_raw_payload='x'*1000000)}))
            image=supp/'asm.s_150.chimera_evidence.png';image.write_bytes(b'fixture')
            generate('asm',calls,ctx,root/'out',supplement=supp)
            report=(root/'out/report.md').read_text(encoding='utf-8')
            self.assertLess(len(report),15000)
            self.assertNotIn('large_raw_payload',report)
            self.assertIn('scaffold_evidence/'+image.name,report)
            self.assertIn('Internal to one native primary contig',report)
            with (root/'out/review.tsv').open(encoding='utf-8') as handle:row=next(csv.DictReader(handle,delimiter='\t'))
            self.assertIn('chr7 → chr12',row['evidence_summary'])
            self.assertEqual(row['review_disposition'],'PENDING')

    def test_retained_measurement_has_no_cut_proposal(self):
        from chimera_markdown import assess_transitions
        row=dict(id='C01',source_candidate_id='a',scaffold='s',review_priority='Evidence favors retaining this sampled boundary',chromosome_context='chr12 → chr12',cut_bp=10)
        assess_transitions([row],{'a':dict(verified_gap=True)}, {}, {'a':dict(left_chrom='chr4',right_chrom='chr12')})
        self.assertEqual(row['action'],'RETAIN');self.assertEqual(row['cut_bp'],'')
        self.assertEqual(families([row]),[])


if __name__=='__main__':unittest.main()
