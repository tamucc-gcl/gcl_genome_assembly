"""Pinned one-line regional grouping correction and source-coordinate integrity checks."""
import argparse
from collections import defaultdict
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import xml.etree.ElementTree as ET

from pangenome_source_audit import read_table, write_table, union_bp, low_retention
from pangenome_assignment_replay import (identity_aliases, split_command, summarize_split,
    fragment, parse_decisions)

COMMIT='f9b999e54c49a05ad6e9d44e1832c4a8421f28f2'
PATCH='regional-reference-grouping-v1'
OLD='if (interval.value.first != ref || clip_candidates.empty()) {'
NEW='if (interval.value.second != ref || clip_candidates.empty()) {'
CASE='regional_1Mb_fixed'


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def verified_source(source):
    source=Path(source);manifest=json.loads((source/'manifest.json').read_text(encoding='utf-8'))
    if manifest['commit']!=COMMIT:raise ValueError('Unexpected regional tool source commit')
    for name,digest in manifest['files'].items():
        if Path(name).name!=name or sha(source/name)!=digest:raise ValueError('Pinned source checksum mismatch: '+name)
    text=(source/'rgfa-split.cpp').read_text(encoding='utf-8')
    if text.count(OLD)!=1:raise ValueError('Regional patch context must occur exactly once')
    return manifest,text


def cli_regions(directory,lengths):
    result=[]
    for p in Path(directory).glob('split_*.fa_contigs'):
        chrom=p.name[len('split_'):-len('.fa_contigs')]
        for name in p.read_text(encoding='utf-8').splitlines():
            m=re.fullmatch(r'(.*):(\d+)-(\d+)',name)
            q,lo,hi=(m[1],int(m[2])-1,int(m[3])) if m else (name,0,lengths.get(name,0))
            # The CLI also lists the reference chromosome itself in its bin.
            if q == chrom and q not in lengths:continue
            if q not in lengths:raise ValueError('Unexpected fixture query')
            result.append(dict(assembly='fixture',sequence=q,chromosome=chrom,start_0based=lo,end_exclusive=hi))
    return result


def coordinate_audit(regions,sources,ambiguous='_AMBIGUOUS_'):
    lengths={(r['assembly'],r['sequence']):int(r['source_bp']) for r in sources}
    if len(lengths)!=len(sources):raise ValueError('Duplicate source identity')
    by=defaultdict(list)
    for r in regions:
        key=r['assembly'],r['sequence'];lo,hi=int(r['start_0based']),int(r['end_exclusive'])
        if key not in lengths or not 0<=lo<hi<=lengths[key]:raise ValueError('Assigned region outside original source')
        by[key].append((lo,hi,r['chromosome']))
    rows=[]
    for key,length in sorted(lengths.items()):
        spans=sorted(by[key]);end=-1
        for lo,hi,chrom in spans:
            if lo<end:raise ValueError('Overlapping or duplicate source assignment: '+str(key))
            end=hi
        assigned=[(lo,hi) for lo,hi,c in spans if c!=ambiguous]
        amb=[(lo,hi) for lo,hi,c in spans if c==ambiguous]
        rows.append(dict(assembly=key[0],sequence=key[1],source_bp=length,
            assigned_union_bp=union_bp(assigned),assigned_sum_bp=sum(hi-lo for lo,hi in assigned),
            ambiguous_bp=union_bp(amb),unreported_bp=length-union_bp([(lo,hi) for lo,hi,c in spans]),
            duplicate_or_overlapping_bp=0,status='PASS'))
    return rows


def native_regressions(original,fixed,root):
    """Real CLI regression independent of Cactus/tool doubles and biological names."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    results=[]
    for mirror in (False,True):
        q='source_reverse' if mirror else 'source_forward';length=8000000
        blocks=[(0,600000,'refA'),(700000,1300000,'refA'),(1400000,3400000,'refA'),(3500000,8000000,'refB')]
        if mirror:blocks=sorted((length-hi,length-lo,ref) for lo,hi,ref in blocks)
        paf=root/(q+'.paf')
        with paf.open('w',encoding='utf-8') as out:
            for lo,hi,ref in blocks:
                n=hi-lo
                out.write(f'{q}\t{length}\t{lo}\t{hi}\t+\t{ref}\t10000000\t0\t{n}\t{n}\t{n}\t60\n')
            out.write('ordinary\t4500000\t0\t4500000\t+\trefB\t10000000\t0\t4500000\t4500000\t4500000\t60\n')
        coverages={}
        for label,binary in [('original',original),('fixed',fixed)]:
            directory=root/(q+'_'+label);directory.mkdir()
            cmd=[str(Path(binary).resolve()),'-p',str(paf.resolve()),'-b',str((directory/'split_').resolve()),
                 '-n','0.25','-Q','2','-u','1000000','-s','-a','_AMBIGUOUS_','-L',str((directory/'tool.log').resolve())]
            with (directory/'stderr.txt').open('w') as err:subprocess.run(cmd,check=True,stdout=err,stderr=err)
            regions=cli_regions(directory,{q:length,'ordinary':4500000})
            coordinate_audit(regions,[dict(assembly='fixture',sequence=q,source_bp=length),
                                      dict(assembly='fixture',sequence='ordinary',source_bp=4500000)])
            coverages[label]=sum(int(r['end_exclusive'])-int(r['start_0based']) for r in regions
                if r['sequence']==q and r['chromosome']=='refA')
            ordinary=sum(int(r['end_exclusive'])-int(r['start_0based']) for r in regions
                if r['sequence']=='ordinary' and r['chromosome']=='refB')
            if ordinary!=4500000:raise ValueError('Ordinary chromosome native regression failed')
        if coverages['fixed']<=coverages['original'] or coverages['fixed']<3400000:
            raise ValueError('Real CLI regression did not demonstrate regional grouping recovery')
        results.append(dict(fixture=q,**coverages,status='PASS'))
    return results


def build(source,out):
    source=Path(source).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    manifest,text=verified_source(source)
    original=out/'original';corrected=out/'corrected'
    for dest in (original,corrected):
        dest.mkdir()
        for name in manifest['files']:shutil.copyfile(source/name,dest/name)
    (corrected/'rgfa-split.cpp').write_text(text.replace(OLD,NEW),encoding='utf-8')
    patch=''.join(difflib.unified_diff(text.splitlines(True),text.replace(OLD,NEW).splitlines(True),
        fromfile='upstream/rgfa-split.cpp',tofile='corrected/rgfa-split.cpp'))
    (out/'regional-grouping.diff').write_text(patch,encoding='utf-8')
    compiler=shlex.split(os.environ.get('CXX','c++'))
    commands=[]
    for directory in (original,corrected):
        cmd=compiler+shlex.split(os.environ.get('CPPFLAGS',''))+shlex.split(os.environ.get('CXXFLAGS',''))+[
            '-O2','-std=c++14','-include','limits','-include','cstdint',
            '-pthread','-static-libstdc++','-static-libgcc','-I',str(directory),
            str(directory/'rgfa-split_main.cpp'),str(directory/'rgfa-split.cpp'),str(directory/'pafcoverage.cpp'),
            '-o',str(directory/'rgfa-split')]+shlex.split(os.environ.get('LDFLAGS',''))
        commands.append(cmd)
        with (directory/'compile.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
    tests=native_regressions(original/'rgfa-split',corrected/'rgfa-split',out/'regressions')
    shutil.copyfile(corrected/'rgfa-split',out/'rgfa-split-fixed')
    shutil.copymode(corrected/'rgfa-split',out/'rgfa-split-fixed')
    version=subprocess.check_output(compiler+['--version'],text=True).splitlines()[0]
    data=dict(patch=PATCH,upstream_commit=COMMIT,source_manifest=manifest,
        patched_source_sha256=sha(corrected/'rgfa-split.cpp'),patched_binary_sha256=sha(out/'rgfa-split-fixed'),
        original_binary_sha256=sha(original/'rgfa-split'),compiler=version,commands=commands,regressions=tests,
        limits='One grouping comparison changed; assignment thresholds, lengths and uniqueness options unchanged')
    (out/'build.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    write_table(out/'native_regressions.tsv',tests,list(tests[0]))


def audit_paf_membership(directory,regions,selection,identities):
    lengths={(r['assembly'],r['sequence']):int(r['source_bp']) for r in selection}
    aliases=identity_aliases(read_table(identities));allowed=set()
    for r in regions:allowed.add((r['assembly'],r['sequence'],r['chromosome'],int(r['start_0based']),int(r['end_exclusive'])))
    for p in Path(directory).glob('split_*.paf'):
        chrom=p.name[len('split_'):-len('.paf')]
        with p.open(encoding='utf-8') as handle:
            for line in handle:
                f=line.rstrip().split('\t');value=fragment(f[0],aliases,lengths)
                if value is None:
                    if not f[0].startswith(('id=_MINIGRAPH_|','_MINIGRAPH_#')):raise ValueError('Unknown output biological query')
                    continue
                key,offset,end=value
                if (*key,chrom,offset,end) not in allowed:raise ValueError('Mapping does not belong to an assigned source fragment')
                if int(f[1])!=end-offset or not 0<=int(f[2])<=int(f[3])<=end-offset:
                    raise ValueError('Mapping outside assigned source fragment')


def assign(prepared,binary,build_info,contigs,out):
    p=Path(prepared);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    info=json.loads(Path(build_info).read_text(encoding='utf-8'))
    if info['patch']!=PATCH or sha(binary)!=info['patched_binary_sha256']:raise ValueError('Patched binary provenance mismatch')
    provenance=json.loads((p/'preparation.json').read_text(encoding='utf-8'))
    temp=out/'split_work';temp.mkdir()
    cmd=split_command(p/'cactus.sv.gfa',p/'filtered.paf',p/'effective_config.xml',provenance['reference'],
        contigs,'regional_1Mb',temp/'split_',out/'assignment_log.txt')
    cmd[0]=str(Path(binary).resolve())
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n',encoding='utf-8')
    with (out/'tool_stderr.txt').open('w') as err:subprocess.run(cmd,check=True,stdout=err,stderr=err)
    ambiguous=ET.parse(p/'effective_config.xml').getroot().find('graphmap_split').get('ambiguousName','_AMBIGUOUS_')
    rows,regions=summarize_split(temp,p/'selection.tsv',p/'identities.tsv',CASE,ambiguous,out/'surviving_intervals.tsv.gz')
    sources=read_table(p/'selection.tsv')
    audits=coordinate_audit(regions,sources,ambiguous)
    audit_paf_membership(temp,regions,sources,p/'identities.tsv')
    write_table(out/'coordinate_audit.tsv',audits,list(audits[0]))
    write_table(out/'assignments.tsv',rows,list(rows[0]))
    write_table(out/'regions.tsv',regions,['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
    decisions=parse_decisions(out/'assignment_log.txt',identity_aliases(read_table(p/'identities.tsv')))
    write_table(out/'decisions.tsv',decisions,['assembly','sequence','decision','assigned_chromosome','coverage_fraction',
        'min_coverage','uniqueness_ratio','min_uniqueness','reason'])
    shutil.copyfile(build_info,out/'build.json')
    (out/'case.json').write_text(json.dumps(dict(policy=provenance['policy'],case=CASE,coordinate_integrity='PASS'))+'\n',encoding='utf-8')
    shutil.rmtree(temp)


def report(original,fixed,filters,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    old={};new={};audits=[]
    for dirs,target in ((original,old),(fixed,new)):
        policies=set()
        for d in dirs:
            p=Path(d);case=json.loads((p/'case.json').read_text(encoding='utf-8'))
            expected=CASE if target is new else 'regional_1Mb'
            if case['case']!=expected or case['policy'] in policies:raise ValueError('Invalid regional case/policy matrix')
            policies.add(case['policy'])
            for r in read_table(p/'assignments.tsv'):
                k=case['policy'],r['assembly'],r['sequence']
                if k in target:raise ValueError('Duplicate regional comparison')
                target[k]=r
            if target is new:
                if case.get('coordinate_integrity')!='PASS':raise ValueError('Patched coordinate audit missing')
                for r in read_table(p/'coordinate_audit.tsv'):audits.append(dict(policy=case['policy'],**r))
        if policies!={'native','chromosome_isolated','disabled'}:raise ValueError('Incomplete regional policy matrix')
    if not old or old.keys()!=new.keys():raise ValueError('Incomplete regional comparison matrix')
    audit_keys=[(r['policy'],r['assembly'],r['sequence']) for r in audits]
    if len(set(audit_keys))!=len(audit_keys) or set(audit_keys)!=set(new):raise ValueError('Incomplete coordinate audit matrix')
    for r in audits:
        assigned=int(r['assigned_union_bp'])
        if (r['status']!='PASS' or int(r['duplicate_or_overlapping_bp'])!=0 or
            assigned!=int(r['assigned_sum_bp']) or min(assigned,int(r['ambiguous_bp']),int(r['unreported_bp']))<0 or
            assigned+int(r['ambiguous_bp'])+int(r['unreported_bp'])!=int(r['source_bp'])):
            raise ValueError('Invalid coordinate audit accounting')
    selections=read_table(Path(filters[0])/'selection.tsv');rows=[]
    for key,r in sorted(new.items()):
        before=old[key];rows.append(dict(policy=key[0],assembly=key[1],sequence=key[2],
            input_filtered_bp=r['filtered_query_covered_bp'],native_regional_bp=before['surviving_query_covered_bp'],
            fixed_regional_bp=r['surviving_query_covered_bp'],
            delta_bp=int(r['surviving_query_covered_bp'])-int(before['surviving_query_covered_bp']),
            native_chromosomes=before['assigned_chromosomes'],fixed_chromosomes=r['assigned_chromosomes'],
            ambiguous_bp=r['ambiguous_span_bp']))
    write_table(out/'regional_assignment_comparison.tsv',rows,list(rows[0]))
    write_table(out/'regional_coordinate_audit.tsv',audits,list(audits[0]))
    by={(r['policy'],r['assembly'],r['sequence']):r for r in rows}
    controls={(r['assembly'],r['sequence']) for r in selections if r['chromosome_scale']=='yes'
        and r['status']=='represented_in_clip' and '+' not in r['sequence'] and not low_retention(r)}
    decreases=[r for r in rows if (r['assembly'],r['sequence']) in controls and int(r['delta_bp'])<0]
    write_table(out/'regional_control_decreases.tsv',decreases,list(rows[0]))
    text=['# Corrected regional chromosome assignment','','The only upstream source change compares reference chromosome IDs when grouping adjacent evidence intervals. '
        'Coverage, minimum chunk and uniqueness options are unchanged. Production graph construction has not run.','',
        '| Source | Filtering | Input mapped, Mb | Native regional, Mb | Corrected regional, Mb | Corrected chromosomes | Ambiguous, Mb |',
        '|---|---|---:|---:|---:|---|---:|']
    for s in selections:
        if not (low_retention(s) or ('+' in s['sequence'] and int(s['source_bp'])>=10000000)):continue
        for policy in ('native','chromosome_isolated','disabled'):
            r=by[policy,s['assembly'],s['sequence']]
            text.append(f"| {s['assembly']} / {s['sequence']} | {policy} | {int(r['input_filtered_bp'])/1e6:.3f} | "
                f"{int(r['native_regional_bp'])/1e6:.3f} | {int(r['fixed_regional_bp'])/1e6:.3f} | "
                f"{r['fixed_chromosomes'] or 'None'} | {int(r['ambiguous_bp'])/1e6:.3f} |")
    text += ['', '## Source-coordinate integrity', '',
        f'{len(audits)} source/policy audits passed. Assigned regions, including ambiguous regions, are in bounds and non-overlapping. '
        'Every output mapping was checked against its assigned fragment before temporary PAFs were removed. '
        'Mapped coverage uses source-coordinate unions. Unreported and ambiguous bases remain explicit; they are not duplicated into another chromosome.', '',
        f'{len(controls)} ordinary chromosome controls; {len(decreases)} control/policy combinations have less mapped coverage after correction.', '',
        '[Every source comparison](regional_assignment_comparison.tsv) · [Coordinate audit](regional_coordinate_audit.tsv) · '
        '[Control decreases](regional_control_decreases.tsv)', '',
        'Per-assignment folders retain chromosome-region coordinates, mapping intervals, native logs and the exact build provenance. '
        'These checks certify intermediate source-coordinate integrity, not the final graph. A production rebuild still needs '
        'FULL/CLIP source representation and distinct-node accounting audits before final acceptance.','']
    (out/'regional_assignment_report.md').write_text('\n'.join(text),encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True)
    b=s.add_parser('build');b.add_argument('--source',required=True);b.add_argument('--out',required=True)
    a=s.add_parser('assign')
    for k in ('prepared','binary','build_info','out'):a.add_argument('--'+k.replace('_','-'),required=True)
    a.add_argument('--contigs',nargs='+',required=True)
    r=s.add_parser('report')
    for k in ('original','fixed','filters'):r.add_argument('--'+k,nargs='+',required=True)
    r.add_argument('--out',required=True)
    args=vars(p.parse_args());mode=args.pop('mode');{'build':build,'assign':assign,'report':report}[mode](**args)


if __name__=='__main__':main()
