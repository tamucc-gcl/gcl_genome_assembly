#!/usr/bin/env python3
"""Create real FASTA/BAM/PAF/pairs inputs for an isolated Nextflow smoke test.

No biological claim is validated by these synthetic continuous junctions.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import random
import shlex
import subprocess

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('join_investigation',ROOT/'scripts/comparisons/join_investigation.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def main():
    import numpy as np
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out=a.out.resolve();a.out.mkdir(exist_ok=False)
    rng=random.Random(517);size=1_200_000;sequence=''.join(rng.choices('ACGT',k=size));source=a.out/'source.fa'
    source.write_text('>component\n'+sequence+'\n');assembly=a.out/'assessment.fa';assembly.write_text('>s\n'+sequence+'\n')
    agp=a.out/'assembly.agp';agp.write_text(f's\t1\t{size}\t1\tW\tcomponent\t1\t{size}\t+\n')
    libraries=a.out/'libraries.tsv';libraries.write_text('library_id\tqname_prefix\nsingle\tL_\n')
    pairs=a.out/'contacts.pairs';matrix=np.zeros((12,12),dtype=np.uint32);rows=[]
    # Both within-bin and inter-bin pairs exercise exact aggregation/replay.
    for i in range(12):
        for j in range(i,min(12,i+5)):
            for _ in range(200):
                x=i*100000+rng.randrange(100000);y=j*100000+rng.randrange(100000)
                rows.append(f'L_{len(rows)}\tcomponent\t{x+1}\tcomponent\t{y+1}\t+\t-\tUU\n')
                matrix[i,j]+=1
                if i!=j:matrix[j,i]+=1
    # Observable 5 kb bins exercise the production unique-arm contact baseline.
    for start in range(0,size,5000):
        for k in range(10):
            x=start+100+k;y=start+200+k
            rows.append(f'L_{len(rows)}\tcomponent\t{x+1}\tcomponent\t{y+1}\t+\t-\tUU\n')
            matrix[start//100000,start//100000]+=1
    pairs.write_text('#pairs format v1.0\n'+''.join(rows))
    contacts=[]
    for i in range(3,10):contacts.append(dict(position=i*100000,**m.window_counts(matrix,i,3)))
    sam=a.out/'reads.sam'
    with sam.open('w') as dst:
        dst.write(f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:s\tLN:{size}\n')
        for i in range(12):
            start=580000+i*100;seq=sequence[start:start+40000]
            dst.write(f'read{i}\t0\ts\t{start+1}\t60\t40000=\t*\t0\t0\t{seq}\t*\tNM:i:0\n')
    bam=a.out/'hifi.bam';subprocess.run(['samtools','view','-b','-o',str(bam),str(sam)],check=True);subprocess.run(['samtools','index',str(bam)],check=True)
    source_packet=a.out/'source_packet';source_packet.mkdir();final=a.out/'final';final.mkdir();focus=[];mapping_manifest={}
    catalog={peer:dict(id=peer,sample=peer,eligible=True,taxid='1',chromosome_labels={'x':'chr1','y':'chr2'}) for peer in ('B','C')}
    catalog_file=a.out/'catalog.json';m.save(catalog_file,catalog)
    for focal in ('A_hap1','A_hap2'):
        focal_bam=bam
        focal_pairs=pairs;focal_count=len(rows);focal_contacts=contacts
        if focal.endswith('hap2'):
            # A negative assay: reads on both sides, none bridging the junction.
            negative=a.out/'nonbridging.sam'
            with negative.open('w') as dst:
                dst.write(f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:s\tLN:{size}\n')
                for i in range(12):
                    start=(540000 if i<6 else 610000)+(i%6)*100
                    seq=sequence[start:start+40000]
                    dst.write(f'nonbridge{i}\t0\ts\t{start+1}\t60\t40000=\t*\t0\t0\t{seq}\t*\tNM:i:0\n')
            focal_bam=a.out/'nonbridging.bam'
            subprocess.run(['samtools','view','-b','-o',str(focal_bam),str(negative)],check=True)
            subprocess.run(['samtools','index',str(focal_bam)],check=True)
            negative_rows=[r for r in rows if not (int(r.split('\t')[2])<=600000<int(r.split('\t')[4]))]
            focal_pairs=a.out/'nonbridging.pairs';focal_pairs.write_text('#pairs format v1.0\n'+''.join(negative_rows))
            focal_count=len(negative_rows);negative_matrix=matrix.copy();negative_matrix[:6,6:]=0;negative_matrix[6:,:6]=0
            focal_contacts=[dict(position=i*100000,**m.window_counts(negative_matrix,i,3)) for i in range(3,10)]
        paths=[]
        for peer in ('B','C'):
            paf=a.out/(focal+'_'+peer+'.paf');paf.write_text(
                's\t1200000\t0\t600000\t+\tx\t600000\t0\t600000\t600000\t600000\t60\tcg:Z:600000=\n'
                's\t1200000\t600000\t1200000\t+\ty\t600000\t0\t600000\t600000\t600000\t60\tcg:Z:600000=\n')
            paths.append(dict(a=focal,b=peer,path=str(paf)))
        discovery_task=a.out/(focal+'.discover_task');discovery=discovery_task/(focal+'.discovery');discovery.mkdir(parents=True)
        discovery_file=discovery/'discovery.json';m.save(discovery_file,{})
        staged_alignments=[dict(r,path=f'alignments{i:02d}/{Path(r["path"]).name}') for i,r in enumerate(paths)]
        (discovery_task/'.command.sh').write_text('python3 misassembly_discover.py discover --catalog catalog.json --alignments '+shlex.quote(json.dumps(staged_alignments)))
        (discovery_task/'.command.run').write_text('nxf_stage() {\n'+
            'cp -fRL '+shlex.quote(str(catalog_file))+' catalog.json\n'+''.join(
            'mkdir -p '+shlex.quote(str(Path(s['path']).parent))+' && cp -fRL '+shlex.quote(r['path'])+' '+shlex.quote(s['path'])+'\n'
            for r,s in zip(paths,staged_alignments))+'}\n')
        mapping=a.out/(focal+'.mapping.json');m.save(mapping,dict(assessment_sha256=m.sha(assembly),assembly=focal,sample='A',coordinate_stage='pre_finishing'))
        assessment_task=a.out/(focal+'.assessment_task');retained=assessment_task/(focal+'.evidence');retained.mkdir(parents=True)
        args=dict(fasta=assembly,discovery=discovery_file,mapping=mapping,agp=agp,pairs=focal_pairs,source=source,libraries=libraries)
        staged={key:f'{key}/{Path(value).name}' for key,value in args.items()}
        (assessment_task/'.command.sh').write_text('python3 misassembly_assess.py '+' '.join('--'+k+' '+shlex.quote(v) for k,v in staged.items()))
        # Match Crest stageInMode=copy: staged copies disappear with scratch,
        # while the retained wrapper records their exact upstream sources.
        (assessment_task/'.command.run').write_text('nxf_stage() {\n'+''.join(
            'mkdir -p '+shlex.quote(str(Path(v).parent))+' && cp -fRL '+shlex.quote(str(args[k]))+' '+shlex.quote(v)+'\n' for k,v in staged.items())+'}\n')
        event=dict(id='C07',scaffold='s',left='chr1',right='chr2',lo=599000,hi=601000,
                   alternatives=[dict(id='C07R01',lo=584000,hi=586000,probe_bp=585000)],cut_options=[],contact_profiles={'single':focal_contacts})
        data=dict(assembly=focal,sample='A',coordinate_stage='pre_finishing',assessment_fasta=str(assembly),assessment_sha256=m.sha(assembly),events=[event],
                  alignments=[dict(a=r['a'],b=r['b'],sha256=m.sha(r['path'])) for r in paths],
                  hic_audit=dict(source_sha256=m.sha(source),agp_sha256=m.sha(agp),audit=dict(input_pairs=focal_count)))
        m.save(retained/'evidence.json',data);pub=source_packet/'assemblies'/focal;pub.mkdir(parents=True);m.save(pub/'evidence.json',data)
        mapping_manifest[focal]=dict(bam=str(focal_bam),index=str(focal_bam)+'.bai',provenance=str(retained/'evidence.json'))
        (final/(focal+'.name_map.tsv')).write_text('old_name\tnew_name\tlength\tclass\ns\tchr1_1+chr2_1\t1200000\tcomposite\n')
        focus.append(dict(assembly=focal,event='C07',role='target'))
    m.save(source_packet/'mapping_manifest.json',mapping_manifest);m.write(a.out/'focus.tsv',focus,['assembly','event','role'])
    print(a.out)


if __name__=='__main__':main()
