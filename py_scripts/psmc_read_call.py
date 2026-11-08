"""Call complete scaffolds concurrently, preserving the serial calling options."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import time


def selected_scaffolds(fai,minimum):
    rows=[]
    for line in Path(fai).read_text().splitlines():
        name,length,*_=line.split('\t')
        if int(length)>=minimum:rows.append((name,int(length)))
    if not rows:raise ValueError('No reference scaffolds meet the PSMC minimum length')
    return rows


def run(a):
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    rows=selected_scaffolds(str(a.reference)+'.fai',a.min_contig_bp)
    if min(a.cpus,a.workers,a.max_pileup_depth)<1:raise ValueError('Positive CPUs, workers and maximum pileup depth required')
    workers=min(len(rows),a.workers,max(1,a.cpus//2))
    shards=out/'shards';shards.mkdir(exist_ok=True)
    with (out/'depth.tsv').open('w') as f:
        subprocess.run(['samtools','coverage','-q',str(a.mapq),'-Q',str(a.baseq),a.bam],stdout=f,check=True)
    def call(item):
        index,(name,length)=item;start=time.monotonic();stem=shards/f'{index:06d}'
        bed=stem.with_suffix('.bed');bed.write_text(f'{name}\t0\t{length}\n')
        vcf=stem.with_suffix('.vcf.gz');log=stem.with_suffix('.log')
        pileup=['bcftools','mpileup','-Ou','-f',a.reference,'-q',str(a.mapq),'-Q',str(a.baseq),
            '-R',str(bed),'-d',str(a.max_pileup_depth),'-a','FORMAT/DP,FORMAT/AD',a.bam]
        caller=['bcftools','call','--ploidy','2','-m','-A','--keep-unseen-allele','-a','GQ','-Oz','-o',str(vcf)]
        with log.open('w') as f:
            producer=subprocess.Popen(pileup,stdout=subprocess.PIPE,stderr=f)
            try:
                result=subprocess.run(caller,stdin=producer.stdout,stderr=f)
                producer.stdout.close()
                if result.returncode:
                    producer.terminate();producer.wait()
                    raise subprocess.CalledProcessError(result.returncode,caller)
                code=producer.wait()
                if code:raise subprocess.CalledProcessError(code,pileup)
            finally:
                if producer.poll() is None:producer.kill();producer.wait()
        print(f'Completed scaffold {index+1}/{len(rows)}: {name}',flush=True)
        return dict(index=index,scaffold=name,length_bp=length,seconds=time.monotonic()-start,
                    pileup=pileup,call=caller,vcf=str(vcf),log=str(log))
    # Whole scaffolds avoid artificial within-scaffold calling boundaries. Sort
    # by FAI index before concatenation, regardless of worker completion order.
    with ThreadPoolExecutor(max_workers=workers) as pool:completed=list(pool.map(call,enumerate(rows)))
    listing=shards/'vcfs.txt';listing.write_text(''.join(r['vcf']+'\n' for r in completed))
    subprocess.run(['bcftools','concat','-f',str(listing),'-Oz','-o',str(out/'diploid.vcf.gz')],check=True)
    subprocess.run(['bcftools','index','--tbi',str(out/'diploid.vcf.gz')],check=True)
    with (out/'calling.log').open('w') as dest:
        for r in completed:dest.write('## '+r['scaffold']+'\n'+Path(r['log']).read_text()+'\n')
    (out/'calling_manifest.json').write_text(json.dumps(dict(workers=workers,requested_workers=a.workers,
        cpus=a.cpus,method='Whole-scaffold parallel pileup and diploid calling; FAI-order concatenation',
        filters=dict(mapq=a.mapq,baseq=a.baseq,max_pileup_depth=a.max_pileup_depth,min_contig_bp=a.min_contig_bp),
        scaffolds=completed),indent=2)+'\n')
    (out/'NO_CONFIDENT_BED').touch()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['reference','bam','out']:p.add_argument('--'+name,required=True)
    for name in ['min-contig-bp','mapq','baseq','max-pileup-depth','cpus','workers']:p.add_argument('--'+name,type=int,required=True)
    run(p.parse_args())


if __name__=='__main__':main()
