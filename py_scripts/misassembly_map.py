#!/usr/bin/env python3
"""Map once per assessment assembly or validate an explicitly retained BAM."""
import argparse
import json
import subprocess
from pathlib import Path
from misassembly_core import sha


def main():
    p=argparse.ArgumentParser()
    for name in ('fasta','discovery','out'):p.add_argument('--'+name,required=True)
    sources=p.add_mutually_exclusive_group();sources.add_argument('--reads');sources.add_argument('--bam')
    p.add_argument('--index');p.add_argument('--provenance');p.add_argument('--threads',type=int,default=8)
    a=p.parse_args();out=Path(a.out);out.mkdir();discovery=json.loads(Path(a.discovery).read_text())
    checksum=sha(a.fasta)
    if discovery['assessment_sha256']!=checksum:raise ValueError('Discovery/FASTA checksum mismatch')
    record=dict(assessment_sha256=checksum,coordinate_stage='pre_finishing',assembly=discovery['assembly'],sample=discovery['sample'],status='not_needed')
    if discovery['events'] and (a.reads or a.bam):
        fasta=out/'assessment.fa';fasta.symlink_to(Path(a.fasta).resolve())
        subprocess.run(['samtools','faidx',str(fasta)],check=True)
        sizes={f[0]:int(f[1]) for f in (l.split('\t') for l in Path(str(fasta)+'.fai').read_text().splitlines())}
        if sum(sizes.values())>=8_000_000_000:raise ValueError('Assembly exceeds single minimap2 index budget')
        bam=out/'hifi.bam'
        if a.bam:
            if not a.provenance:raise ValueError('Retained BAM requires its evidence provenance')
            provenance=json.loads(Path(a.provenance).read_text())
            if provenance.get('assessment_sha256')!=checksum or provenance.get('coordinate_stage')!='pre_finishing' or provenance.get('assembly')!=discovery['assembly'] or provenance.get('sample')!=discovery['sample']:raise ValueError('Reused BAM has different assembly/sample provenance or coordinates')
            bam.symlink_to(Path(a.bam).resolve())
            index=Path(a.index) if a.index else Path(str(a.bam)+'.bai')
            if not index.exists():raise ValueError('Retained BAM requires adjacent .bai')
            Path(str(bam)+'.bai').symlink_to(index.resolve())
            record['status']='reused'
        else:
            if not a.reads:raise ValueError('Missing HiFi reads')
            # Stream to sorting: do not create a giant unused SAM output.
            command=['minimap2','-ax','map-hifi','-I','8G','--secondary=yes','-N','50','-p','0.5','-t',str(a.threads),str(fasta),a.reads]
            record['mapping_command']=command;record['reads_bytes']=Path(a.reads).stat().st_size
            mapper=subprocess.Popen(command,stdout=subprocess.PIPE)
            try:
                subprocess.run(['samtools','sort','-@','2','-m','2G','-o',str(bam),'-'],stdin=mapper.stdout,check=True)
            except BaseException:
                mapper.terminate();raise
            finally:
                mapper.stdout.close();status=mapper.wait()
            if status:raise RuntimeError('HiFi mapping failed')
            subprocess.run(['samtools','index',str(bam)],check=True);record['status']='mapped'
        subprocess.run(['samtools','quickcheck',str(bam)],check=True)
        idx=subprocess.check_output(['samtools','idxstats',str(bam)],text=True)
        if {f[0]:int(f[1]) for f in (l.split('\t') for l in idx.splitlines()) if f[0]!='*'}!=sizes:raise ValueError('BAM dictionary differs from assessment FASTA')
        # Nextflow moves task outputs off node-local scratch. Paths inside this
        # packet must be relative until the consumer resolves its staged packet.
        record.update(bam='hifi.bam',index='hifi.bam.bai')
        fasta.unlink();Path(str(fasta)+'.fai').unlink()
    elif discovery['events']:record['status']='unavailable'
    (out/'mapping.json').write_text(json.dumps(record,indent=2),encoding='utf-8')


if __name__=='__main__':main()
