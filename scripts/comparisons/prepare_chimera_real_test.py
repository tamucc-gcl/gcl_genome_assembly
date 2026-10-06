#!/usr/bin/env python3
"""Build explicit current-input manifests for an isolated CHIMERA workflow test."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


def digest(path):
    with open(path,'rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def require(path):
    path=Path(path).resolve()
    if not path.is_file():raise ValueError('Missing required input: '+str(path))
    return str(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assessment',required=True)
    p.add_argument('--bam-root',required=True)
    p.add_argument('--species',required=True)
    p.add_argument('--sample',required=True)
    p.add_argument('--reviewed-template',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    root=Path(a.assessment).resolve();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    harm=root/'assembly/harmonization'
    quality=require(harm/(a.species+'.chromosome_sets.tsv'))
    with open(quality) as handle:quality_rows=list(csv.DictReader(handle,delimiter='\t'))
    reference=require(harm/(a.species+'.reference_id.txt'))
    reference_id=Path(reference).read_text().strip()
    cohort=[];targets=[];bams={}
    for row in quality_rows:
        assembly=row['id'];sample=re.sub(r'_hap\d+$','',assembly)
        fasta=require(root/'assembly/scaffold/yahs_round2'/(assembly+'_round2_scaffolds.fa'))
        names=require(harm/(assembly+'.harmonized_name_map.tsv'))
        record=dict(meta=dict(id=assembly,sample=sample,taxid=a.species,assembler='hifiasm',hic=True),fasta=fasta,names=names)
        cohort.append(record)
        if sample==a.sample:
            record.update(agp=require(root/'assembly/scaffold/yahs_round2'/(assembly+'_round2_scaffolds_final.agp')),
                          graph=require(root/'assembly/contig/hifiasm'/(sample+'.'+assembly[len(sample)+1:]+'.p_ctg.gfa')),
                          paf=require(harm/(assembly+'.ref.paf.gz')),
                          pairs=require(root/'bam/hic/scaffold/filtered'/(assembly+'.pairs.gz')),
                          source=require(root/'assembly/scaffold/misassembly_correction'/(assembly+'_corrected.fasta')),
                          libraries=require(root/'bam/hic/scaffold/raw'/(assembly+'.readsets.tsv')))
            folder=Path(a.bam_root)/(assembly+'.sequence_context')
            bam=require(folder/'hifi.bam')
            index=next((str(x.resolve()) for x in [Path(bam+'.bai'),Path(bam+'.csi')] if x.is_file()),None)
            if not index:raise ValueError('Missing retained BAM index: '+bam)
            provenance=require(folder/'provenance.json')
            stamp=json.loads(Path(provenance).read_text())
            if stamp.get('sha256')!=digest(fasta) or stamp.get('coordinate_stage')!='pre_finishing':
                raise ValueError('Retained BAM does not match current assessment: '+assembly)
            bams[assembly]=dict(bam=bam,index=index,provenance=provenance)
            targets.append(record)
    if len(targets)!=2:raise ValueError('This diploid test requires exactly two current target haplotypes')
    if reference_id not in {r['meta']['id'] for r in cohort}:raise ValueError('Reference absent from cohort')
    manifest=dict(species=a.species,reference_id=reference_id,quality=quality,cohort=cohort,targets=targets,
                  candidates=require(harm/(a.species+'.chimera_candidates.tsv')),
                  reference_names=require(harm/(reference_id+'.harmonized_name_map.tsv')))
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'bams.json').write_text(json.dumps(bams,indent=2)+'\n')
    # The reviewed regression action is test data, never discovery logic.
    actions=json.loads(Path(a.reviewed_template).read_text())
    target_by_id={r['meta']['id']:r for r in targets}
    for action in actions:
        record=target_by_id.get(action['assembly'])
        if not record or action['assessment_sha256']!=digest(record['fasta']):
            raise ValueError('Reviewed regression action does not match assessed assembly')
    with (out/'reviewed-actions.tsv').open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=list(actions[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(actions)
    (out/'reviewed-actions.json').write_text(json.dumps(actions,indent=2)+'\n')
    print(json.dumps(dict(status='READY',cohort=len(cohort),targets=len(targets),reference=reference_id)))


if __name__=='__main__':main()
