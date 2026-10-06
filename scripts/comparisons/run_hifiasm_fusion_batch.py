#!/usr/bin/env python3
"""Run one isolated treatment in a parallel hifiasm fusion experiment; never cut assemblies."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'py_scripts'))
from chimera_origin import paf_rows, write_table
from junction_focus import union
from trace_chimera_origins import Fasta
from chimera_origin import table

TREATMENTS = [('both', ('A','B'), 1), ('A_only', ('A',), 1),
              ('B_only', ('B',), 1), ('hifi_only', (), 1), ('both_no_postjoin', ('A','B'), 0)]


def command(sample, threads, hifi, libraries, postjoin):
    # Exact graph-generation/partition settings from CTlk job 1499288.
    args = ['hifiasm', '-o', sample, '-t', str(threads), '-k','51','-w','51','-f','37',
            '-D','5.0','-N','100','-r','3','-z','0','--max-kocc','2000',
            '-a','4','-m','10000000','-p','0','-n','3','-x','0.8','-y','0.2',
            '-u',str(postjoin),'--lowQ','70','--b-cov','0','--h-cov','-1',
            '--m-rate','0.75','--ctg-n','3','-l','3','-s','0.55','-O','1',
            '--n-hap','2','--scaf-gap','3000000','--telo-m','CCCTAA',
            '--telo-p','1','--telo-d','2000','--telo-s','500','--s-base','0.5',
            '--n-weight','3','--n-perturb','10000','--f-perturb','0.1','--l-msjoin','500000']
    if libraries:
        args += ['--h1', ','.join(str(pair[0]) for pair in libraries),
                 '--h2', ','.join(str(pair[1]) for pair in libraries)]
    return args + [str(hifi)]


def execute(args, cwd, logfile):
    print(shlex.join([str(x) for x in args]), flush=True)
    with logfile.open('w') as out:
        subprocess.run(args, cwd=cwd, stdout=out, stderr=subprocess.STDOUT, check=True)


def composition(hits):
    groups = defaultdict(list)
    for h in hits:
        if h['mapq'] >= 20 and h['query_end']-h['query_start'] >= 2000:
            groups[h['query'],h['target']].append(h)
    rows = []
    for (query,target), hs in sorted(groups.items()):
        bases = sum(b-a for a,b in union((h['query_start'],h['query_end']) for h in hs))
        rows.append(dict(contig=query,reference_scaffold=target,contig_length=hs[0]['query_length'],
                         union_query_span_bp=bases,fraction=bases/hs[0]['query_length'],
                         strands=';'.join(sorted({h['strand'] for h in hs})),
                         interpretation='alignment_span_screen_not_unique_chromosome_or_join_validation'))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,required=True)
    p.add_argument('--sample',default='Sde-CTlk_104')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--treatment',type=int,choices=range(5),required=True)
    p.add_argument('--threads',type=int,default=48)
    p.add_argument('--assembler-env',help='Separate existing conda environment for hifiasm/gfatools')
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--reference-name-map',type=Path,required=True)
    p.add_argument('--baseline-hap1',type=Path,required=True)
    p.add_argument('--baseline-hap2',type=Path,required=True)
    p.add_argument('--target-contigs',type=Path,required=True,
                   help='Investigation data: haplotype/contig TSV for targeted original-path comparisons')
    a=p.parse_args()
    if not a.sample or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-' for c in a.sample) or a.sample in ('.','..') or a.threads<=0:
        p.error('Safe sample identifier and positive threads required')
    root=a.input_root.resolve()
    hifi=root/'fastq/hifi'/f'{a.sample}.fastq.gz'
    libs={key: (root/'fastq/hic/trimmed'/f'{a.sample}__{lib}_run1_R1.trim.fastq.gz',
                root/'fastq/hic/trimmed'/f'{a.sample}__{lib}_run1_R2.trim.fastq.gz')
          for key,lib in [('A','Ex2'),('B','Ex3')]}
    reference=a.reference.resolve(); names=a.reference_name_map.resolve()
    baseline=[a.baseline_hap1.resolve(),a.baseline_hap2.resolve()]
    label, selected, postjoin=TREATMENTS[a.treatment]
    inputs=[hifi,reference,names]+baseline+[a.target_contigs.resolve()]+[p for key in selected for p in libs[key]]
    for path in inputs:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f'Missing/nonempty input required: {path}')
    env_prefix=[]
    if a.assembler_env:
        env_prefix=['conda','run','--no-capture-output',
                    '-p' if '/' in a.assembler_env else '-n',a.assembler_env]
    version=subprocess.check_output(env_prefix+['hifiasm','--version'],text=True).strip()
    if version!='0.25.0-r726':
        raise ValueError(f'Experiment requires original hifiasm 0.25.0-r726, found {version}')
    out=a.out.resolve()/label; out.mkdir(parents=True,exist_ok=False)
    (out/'runner.py').write_bytes(Path(__file__).read_bytes())
    manifest=dict(treatment=label,postjoin=postjoin,libraries=list(selected),version=version,
                  inputs=[dict(path=str(p),resolved_path=str(p.resolve()),bytes=p.stat().st_size,
                               mtime_ns=p.stat().st_mtime_ns) for p in inputs],
                  input_identity='Same retained input paths; no repeated whole-read-file hashing in each array task',
                  status='RUNNING',cut_authorized=False)
    (out/'status.json').write_text(json.dumps(manifest,indent=2)+'\n')
    args=env_prefix+command(a.sample,a.threads,hifi,[libs[key] for key in selected],postjoin)
    (out/'command.json').write_text(json.dumps(args,indent=2)+'\n')
    try:
        tool_versions={tool:subprocess.check_output([tool,'--version'],text=True,
                      stderr=subprocess.STDOUT).splitlines()[0] for tool in ('minimap2','samtools')}
        (out/'analysis_versions.json').write_text(json.dumps(tool_versions,indent=2)+'\n')
        queries=[]
        targets=table(a.target_contigs)
        for original_hap, source in enumerate(baseline,1):
            sourcefa=Fasta(source,out/f'original_hap{original_hap}.fa','samtools')
            query=out/f'original_hap{original_hap}.targets.fa'
            with query.open('w') as dest:
                wanted=[r['contig'] for r in targets if int(r['haplotype'])==original_hap]
                if not wanted or len(set(wanted))!=len(wanted):
                    raise ValueError('Target manifest requires distinct contigs for both haplotypes')
                for contig in wanted:
                    dest.write(f'>{contig}\n{sourcefa.fetch(contig,0,sourcefa.lengths[contig])}\n')
            queries.append(query)
        execute(args,out,out/'hifiasm.log')
        prefix='hic' if selected else 'bp'
        (out/'reference_name_map.tsv').write_bytes(names.read_bytes())
        for hap in (1,2):
            gfa=out/f'{a.sample}.{prefix}.hap{hap}.p_ctg.gfa'
            fasta=out/f'hap{hap}.fa'
            with fasta.open('w') as dest:
                subprocess.run(env_prefix+['gfatools','gfa2fa',str(gfa)],stdout=dest,check=True)
            paf=out/f'hap{hap}.reference.paf'
            cmd=['minimap2','-x','asm5','-I','8G','-c','--eqx','--secondary=yes',
                 '-N','50','-p','0.5','-t',str(a.threads),str(reference),str(fasta)]
            with paf.open('w') as dest,(out/f'hap{hap}.reference.log').open('w') as log:
                subprocess.run(cmd,stdout=dest,stderr=log,check=True)
            write_table(out/f'hap{hap}.composition.tsv',
                        ['contig','reference_scaffold','contig_length','union_query_span_bp','fraction','strands','interpretation'],
                        composition(paf_rows(paf)))
            # Compare BOTH original haplotypes to each treatment; hap numbers can swap.
            for original_hap, source in enumerate(queries,1):
                paf=out/f'original_hap{original_hap}.to_hap{hap}.paf'
                cmd=['minimap2','-x','asm5','-I','8G','-c','--eqx','--secondary=yes',
                     '-N','50','-p','0.5','-t',str(a.threads),str(fasta),str(source)]
                with paf.open('w') as dest,(out/f'original_hap{original_hap}.to_hap{hap}.log').open('w') as log:
                    subprocess.run(cmd,stdout=dest,stderr=log,check=True)
        manifest['status']='SUCCESS'
    except Exception as exc:
        manifest.update(status='FAILED',error=str(exc))
        raise
    finally:
        manifest['retained_files']=[dict(name=p.name,bytes=p.stat().st_size)
                                    for p in sorted(out.iterdir()) if p.is_file()]
        (out/'status.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__': main()
