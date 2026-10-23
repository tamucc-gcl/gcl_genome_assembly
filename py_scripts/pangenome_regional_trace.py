"""Observe pinned regional assignment gates; refuse any change to assignment outputs."""
import argparse
from collections import Counter, defaultdict
import csv
import difflib
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import xml.etree.ElementTree as ET

from pangenome_regional_patch import verified_source, OLD, NEW, sha, coordinate_audit, audit_paf_membership
from pangenome_assignment_replay import split_command, summarize_split, identity_aliases
from pangenome_source_audit import read_table, write_table, query_identity, low_retention

TRACE_VERSION='regional-gates-observation-v1'


def once(text,old,new):
    if text.count(old)!=1:raise ValueError('Trace insertion context must occur exactly once: '+old[:80])
    return text.replace(old,new)


def instrument(text):
    """Only insert local statistics and log writes; preserve all original conditions."""
    marker='''        for (auto& query_coverage : coverage_map) {
            vector<CoverageInterval> intervals;'''
    text=once(text,marker,'''        for (auto& query_coverage : coverage_map) {
            for (const auto& trace_ref : query_coverage.second) {
                int64_t trace_max = 0, trace_bp = 0, trace_aligned = 0;
                for (const auto& trace_i : trace_ref.second) {
                    trace_max = std::max(trace_max, trace_i.stop - trace_i.start + 1);
                    trace_bp += trace_i.stop - trace_i.start + 1;
                    trace_aligned += trace_i.value.first;
                }
                log_stream << "TRACE\\tbefore_overlap\\t" << query_coverage.first << "\\t"
                           << contigs[trace_ref.first] << "\\t" << trace_ref.second.size() << "\\t"
                           << trace_max << "\\t" << trace_bp << "\\t" << trace_aligned << endl;
            }
            vector<CoverageInterval> intervals;''')
    marker='''    // total up coverage by reference contig
    map<int64_t, int64_t> coverage_by_contig;'''
    text=once(text,marker,'''    // Observation only: input to smoothing, after overlap resolution.
    map<int64_t, vector<int64_t>> trace_stats;
    for (const auto& trace_i : intervals) {
        auto& trace_s = trace_stats[trace_i.value.second];
        if (trace_s.empty()) trace_s.resize(4, 0);
        ++trace_s[0];
        trace_s[1] = std::max(trace_s[1], trace_i.stop - trace_i.start + 1);
        trace_s[2] += trace_i.stop - trace_i.start + 1;
        trace_s[3] += trace_i.value.first;
    }
    for (const auto& trace_s : trace_stats) {
        log_stream << "TRACE\\tafter_overlap\\t" << query_name << "\\t" << ref_contigs[trace_s.first];
        for (auto trace_v : trace_s.second) log_stream << "\\t" << trace_v;
        log_stream << endl;
    }
    // total up coverage by reference contig
    map<int64_t, int64_t> coverage_by_contig;''')
    marker='''    // look for regions to possibly clip out, by finding runs of consecutive'''
    text=once(text,marker,'''    log_stream << "TRACE\\tdominant\\t" << query_name << "\\t" << ref_contigs[top.first]
               << "\\t" << top.second << "\\t" << min_chunk << endl;
    // look for regions to possibly clip out, by finding runs of consecutive''')
    marker='''        if (max_interval_length > min_chunk &&
            (double)total_coverage / (max_pos - min_pos + 1) >= min_coverage) {'''
    text=once(text,marker,'''        log_stream << "TRACE\\tcandidate\\t" << query_name << "\\t"
                   << ref_contigs[intervals[clip_candidates[i][0]].value.second] << "\\t"
                   << min_pos << "\\t" << (max_pos + 1) << "\\t" << clip_candidates[i].size() << "\\t"
                   << max_interval_length << "\\t" << min_chunk << "\\t" << total_coverage << "\\t"
                   << (double)total_coverage / (max_pos - min_pos + 1) << "\\t"
                   << min_coverage << "\\t" << (max_interval_length > min_chunk) << "\\t"
                   << ((double)total_coverage / (max_pos - min_pos + 1) >= min_coverage) << endl;
        if (max_interval_length > min_chunk &&
            (double)total_coverage / (max_pos - min_pos + 1) >= min_coverage) {''')
    marker='''            log_stream << "Unable to smooth small fragment: "'''
    text=once(text,marker,'''            log_stream << "TRACE\\tsmall_dominant\\t" << query_name << "\\t"
                       << ref_contigs[smooth_intervals[i].value.second] << "\\t" << smooth_intervals[i].start
                       << "\\t" << (smooth_intervals[i].stop + 1) << "\\t" << min_len_filter << endl;
            log_stream << "Unable to smooth small fragment: "''')
    marker='''    if (adjusted_coverage > min_coverage) {'''
    text=once(text,marker,'''    log_stream << "TRACE\\tfinal\\t" << query_name << "\\t" << ref_contigs[top.first]
               << "\\t" << top.second << "\\t" << query_length << "\\t" << masked_bases << "\\t"
               << softclip << "\\t" << total_clip_length << "\\t" << adjusted_coverage << "\\t"
               << min_coverage << "\\t" << (adjusted_coverage > min_coverage) << endl;
    if (adjusted_coverage > min_coverage) {''')
    return text


def build(source,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    manifest,text=verified_source(source);commands=[]
    versions={'native':text,'native_trace':instrument(text),'fixed_trace':instrument(text.replace(OLD,NEW))}
    for name,code in versions.items():
        dest=out/name;dest.mkdir()
        for filename in manifest['files']:shutil.copyfile(source/filename,dest/filename)
        (dest/'rgfa-split.cpp').write_text(code,encoding='utf-8')
        cmd=shlex.split(os.environ.get('CXX','c++'))+shlex.split(os.environ.get('CPPFLAGS',''))+shlex.split(os.environ.get('CXXFLAGS',''))+[
            '-O2','-std=c++14','-include','limits','-include','cstdint','-pthread',
            '-static-libstdc++','-static-libgcc','-I',str(dest),str(dest/'rgfa-split_main.cpp'),
            str(dest/'rgfa-split.cpp'),str(dest/'pafcoverage.cpp'),'-o',str(out/(name+'.bin'))]+shlex.split(os.environ.get('LDFLAGS',''))
        commands.append(cmd)
        with (dest/'compile.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
        if name!='native':
            original=versions['native'] if name=='native_trace' else text.replace(OLD,NEW)
            (out/(name+'.diff')).write_text(''.join(difflib.unified_diff(original.splitlines(True),code.splitlines(True),
                fromfile=name+'/without_observation.cpp',tofile=name+'/with_observation.cpp')),encoding='utf-8')
    regressions=native_gate_regressions(out/'native.bin',out/'native_trace.bin',out/'regressions')
    write_table(out/'gate_regressions.tsv',regressions,list(regressions[0]))
    metadata=dict(version=TRACE_VERSION,source_manifest=manifest,commands=commands,regressions=regressions,
        binaries={n:sha(out/(n+'.bin')) for n in versions},
        scope='Observation only; original decision expressions and thresholds unchanged')
    (out/'trace_build.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')


def output_fingerprints(directory):
    return {p.name:sha(p) for p in sorted(Path(directory).glob('split_*'))
            if p.is_file() and p.suffix in ('.paf','.fa_contigs')}


def verify_equivalent(plain,observed):
    left,right=output_fingerprints(plain),output_fingerprints(observed)
    if not left or left!=right:raise ValueError('Observation changed native assignment output bytes')
    return len(left)


def native_gate_regressions(plain_binary,observed_binary,root):
    """Known gate outcomes from actual upstream CLI, not a Python algorithm copy."""
    root=Path(root);root.mkdir()
    identities=root/'ids.tsv';identities.write_text('id\tgraph_name\nfixture\tTRACE.1\n')
    cases={
        'length_equal':([(0,1000000,'A',1000000),(1000000,8000000,'B',7000000)],'length_failed_runs',1),
        'length_over':([(0,1000001,'A',1000001),(1000001,8000000,'B',6999999)],'locally_accepted_runs',1),
        'low_density':([(0,2000000,'A',200000),(2000000,8000000,'B',6000000)],'density_failed_runs',1),
        'low_whole_query':([(0,2000000,'A',1000000),(2000000,8000000,'B',1100000)],'whole_query_coverage_pass','0'),
        'overlap_removed':([(0,2000000,'A',2000000),(0,8000000,'B',8000000)],'post_overlap_intervals',0)}
    results=[]
    for name,(blocks,field,expected) in cases.items():
        directory=root/name;directory.mkdir();paf=directory/'in.paf'
        paf.write_text(''.join(f'id=TRACE.1|{name}\t8000000\t{lo}\t{hi}\t+\t{chrom}\t10000000\t0\t{hi-lo}\t{matches}\t{hi-lo}\t60\n'
            for lo,hi,chrom,matches in blocks))
        dirs={}
        for label,binary in [('plain',plain_binary),('observed',observed_binary)]:
            dest=directory/label;dest.mkdir();dirs[label]=dest
            cmd=[str(Path(binary).resolve()),'-p',str(paf.resolve()),'-b',str((dest/'split_').resolve()),
                 '-n','.25','-Q','2','-u','1000000','-s','-a','_AMBIGUOUS_','-L',str((dest/'tool.log').resolve())]
            with (dest/'stderr.txt').open('w') as err:subprocess.run(cmd,check=True,stdout=err,stderr=err)
        count=verify_equivalent(dirs['plain'],dirs['observed'])
        rows=parse_trace(dirs['observed']/'tool.log',identities,
            [dict(assembly='fixture',sequence=name,source_bp=8000000)],directory/'events.tsv.gz')
        secondary=next(r for r in rows if r['chromosome']=='A')
        if secondary[field]!=expected:raise ValueError('Native gate regression failed: '+name)
        results.append(dict(fixture=name,measured_field=field,expected=expected,observed=secondary[field],
            compared_files=count,byte_identical='yes',status='PASS'))
    return results


FIELDS={
 'before_overlap':['query','chromosome','intervals','max_interval_bp','interval_span_sum_bp','aligned_bp'],
 'after_overlap':['query','chromosome','intervals','max_interval_bp','interval_span_sum_bp','aligned_bp'],
 'dominant':['query','chromosome','aligned_bp','min_chunk_bp'],
 'candidate':['query','chromosome','start_0based','end_exclusive','intervals','max_interval_bp','min_chunk_bp',
              'aligned_bp','density','min_density','length_pass','density_pass'],
 'small_dominant':['query','chromosome','start_0based','end_exclusive','min_length_bp'],
 'final':['query','chromosome','dominant_aligned_bp','query_bp','masked_bp','softclip_bp','accepted_secondary_span_bp',
          'adjusted_coverage','min_coverage','coverage_pass']}


def parse_trace(log,identities,selection,events_path):
    aliases=identity_aliases(read_table(identities));sources={(r['assembly'],r['sequence']) for r in selection}
    stats=defaultdict(lambda:dict(before={},after={},runs=0,length_fail=0,density_fail=0,
        accepted=0,max_span=0,max_interval=0,min_chunk='',small=0));final={};dominants={}
    header=['stage','assembly','sequence']+list(dict.fromkeys(k for keys in FIELDS.values() for k in keys if k!='query'))
    with open(log,encoding='utf-8') as inp,gzip.open(events_path,'wt',encoding='utf-8',newline='') as dest:
        writer=csv.DictWriter(dest,fieldnames=header,delimiter='\t',lineterminator='\n');writer.writeheader()
        for line in inp:
            if not line.startswith('TRACE\t'):continue
            values=line.rstrip('\n').split('\t');stage=values[1]
            if stage not in FIELDS or len(values[2:])!=len(FIELDS[stage]):raise ValueError('Malformed regional trace event')
            row=dict(zip(FIELDS[stage],values[2:]));key=query_identity(row.pop('query'),aliases)
            if key not in sources:raise ValueError('Unexpected trace query')
            writer.writerow(dict(stage=stage,assembly=key[0],sequence=key[1],**row))
            chrom=row['chromosome'];record=stats[(*key,chrom)]
            if stage in ('before_overlap','after_overlap'):
                field='before' if stage=='before_overlap' else 'after'
                if record[field]:raise ValueError('Duplicate regional trace stage')
                record[field]=row
            elif stage=='candidate':
                record['runs']+=1;record['length_fail']+=row['length_pass']=='0';record['density_fail']+=row['density_pass']=='0'
                record['accepted']+=row['length_pass']=='1' and row['density_pass']=='1'
                record['max_span']=max(record['max_span'],int(row['end_exclusive'])-int(row['start_0based']))
                record['max_interval']=max(record['max_interval'],int(row['max_interval_bp']))
                record['min_chunk']=row['min_chunk_bp']
            elif stage=='small_dominant':record['small']+=1
            elif stage=='dominant':dominants[key]=chrom
            elif stage=='final':
                if key in final:raise ValueError('Duplicate final trace event')
                final[key]=row
    rows=[]
    for key,s in sorted(stats.items()):
        aid,seq,chrom=key;qkey=aid,seq;f=final.get(qkey,{})
        length_fail=s['length_fail'];density_fail=s['density_fail'];accepted=s['accepted']
        dominant=dominants.get(qkey)==chrom
        if not s['after']:reason='All intervals removed by overlap resolution'
        elif not f:reason='Missing final trace'
        elif f['coverage_pass']=='0':reason='Whole-query adjusted coverage failed'
        elif dominant:
            reason='Dominant chromosome; inspect final regions' if not s['small'] else 'Dominant chromosome; small regions removed'
        elif accepted:reason='Secondary region passed local gates; inspect final regions'
        elif s['runs']:
            failures=[]
            if length_fail:failures.append('longest interval too short')
            if density_fail:failures.append('aligned-base density too low')
            reason='Secondary regions rejected: '+'; '.join(failures)
        else:reason='No secondary candidate produced'
        rows.append(dict(assembly=aid,sequence=seq,chromosome=chrom,role='dominant' if dominant else 'secondary',
            input_intervals=s['before'].get('intervals',0),input_max_interval_bp=s['before'].get('max_interval_bp',0),
            input_aligned_bp=s['before'].get('aligned_bp',0),post_overlap_intervals=s['after'].get('intervals',0),
            post_overlap_max_interval_bp=s['after'].get('max_interval_bp',0),candidate_runs=s['runs'],
            largest_candidate_span_bp=s['max_span'],largest_candidate_interval_bp=s['max_interval'],
            min_chunk_bp=s['min_chunk'],length_failed_runs=length_fail,density_failed_runs=density_fail,
            locally_accepted_runs=accepted,removed_small_dominant_regions=s['small'],
            whole_query_coverage_pass=f.get('coverage_pass',''),adjusted_coverage=f.get('adjusted_coverage',''),
            min_coverage=f.get('min_coverage',''),reason=reason))
    if not rows:raise ValueError('No regional trace events')
    if any(r['reason']=='Missing final trace' for r in rows):raise ValueError('Incomplete regional trace')
    return rows


def run(prepared,tools,fixed_binary,fixed_info,contigs,out):
    p=Path(prepared);tools=Path(tools);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    info=json.loads((tools/'trace_build.json').read_text());fixed=json.loads(Path(fixed_info).read_text())
    if info['version']!=TRACE_VERSION or sha(fixed_binary)!=fixed['patched_binary_sha256']:
        raise ValueError('Trace binary provenance mismatch')
    for name,digest in info['binaries'].items():
        if sha(tools/(name+'.bin'))!=digest:raise ValueError('Trace binary checksum mismatch')
    provenance=json.loads((p/'preparation.json').read_text());selection=read_table(p/'selection.tsv')
    ambiguous=ET.parse(p/'effective_config.xml').getroot().find('graphmap_split').get('ambiguousName','_AMBIGUOUS_')
    comparisons=[];all_rows=[];commands=[]
    for variant,plain_binary,observed_binary in [('native',tools/'native.bin',tools/'native_trace.bin'),
                                               ('fixed',Path(fixed_binary),tools/'fixed_trace.bin')]:
        dirs={}
        for label,binary in [('plain',plain_binary),('observed',observed_binary)]:
            directory=out/(variant+'_'+label);directory.mkdir();dirs[label]=directory
            cmd=split_command(p/'cactus.sv.gfa',p/'filtered.paf',p/'effective_config.xml',provenance['reference'],
                contigs,'regional_1Mb',directory/'split_',directory/'tool.log');cmd[0]=str(binary.resolve());commands.append(cmd)
            with (directory/'stderr.txt').open('w') as err:subprocess.run(cmd,check=True,stdout=err,stderr=err)
        files=verify_equivalent(dirs['plain'],dirs['observed'])
        case='trace_'+variant
        assignments,regions=summarize_split(dirs['observed'],p/'selection.tsv',p/'identities.tsv',case,ambiguous)
        audits=coordinate_audit(regions,selection,ambiguous)
        audit_paf_membership(dirs['observed'],regions,selection,p/'identities.tsv')
        write_table(out/(variant+'_coordinate_audit.tsv'),audits,list(audits[0]))
        write_table(out/(variant+'_regions.tsv'),regions,['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
        write_table(out/(variant+'_assignments.tsv'),assignments,list(assignments[0]))
        rows=parse_trace(dirs['observed']/'tool.log',p/'identities.tsv',selection,out/(variant+'_events.tsv.gz'))
        assigned=Counter()
        for r in regions:assigned[r['assembly'],r['sequence'],r['chromosome']]+=int(r['end_exclusive'])-int(r['start_0based'])
        for row in rows:
            row['variant']=variant
            row['final_assigned_span_bp']=assigned[row['assembly'],row['sequence'],row['chromosome']]
            if row['final_assigned_span_bp']>0:row['outcome']='Assigned'
            else:row['outcome']='Not assigned'
        all_rows.extend(rows)
        comparisons.append(dict(variant=variant,compared_output_files=files,byte_identical='yes',coordinate_integrity='PASS'))
        shutil.rmtree(dirs['plain']);shutil.rmtree(dirs['observed'])
    write_table(out/'chromosome_gates.tsv',all_rows,list(all_rows[0]))
    write_table(out/'observation_equivalence.tsv',comparisons,list(comparisons[0]))
    shutil.copyfile(p/'selection.tsv',out/'selection.tsv');shutil.copyfile(tools/'trace_build.json',out/'trace_build.json')
    (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
    (out/'trace_case.json').write_text(json.dumps(dict(policy=provenance['policy'],version=TRACE_VERSION))+'\n')


def report(inputs,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[];checks=[];policies=set();selection=None
    for d in inputs:
        p=Path(d);meta=json.loads((p/'trace_case.json').read_text());policy=meta['policy']
        if policy in policies or meta['version']!=TRACE_VERSION:raise ValueError('Invalid trace policy matrix')
        policies.add(policy)
        for row in read_table(p/'observation_equivalence.tsv'):
            if row['byte_identical']!='yes' or row['coordinate_integrity']!='PASS':raise ValueError('Observation validation failed')
            checks.append(dict(policy=policy,**row))
        if {r['variant'] for r in read_table(p/'observation_equivalence.tsv')}!={'native','fixed'}:raise ValueError('Incomplete observation comparisons')
        rows.extend(dict(policy=policy,**r) for r in read_table(p/'chromosome_gates.tsv'))
        sources=read_table(p/'selection.tsv')
        if selection is None:selection=sources
        elif sorted((r['assembly'],r['sequence'],r['source_bp']) for r in sources)!=sorted(
                (r['assembly'],r['sequence'],r['source_bp']) for r in selection):
            raise ValueError('Inconsistent trace source selections')
    if policies!={'native','chromosome_isolated','disabled'}:raise ValueError('Incomplete trace policies')
    write_table(out/'regional_gate_summary.tsv',rows,list(rows[0]));write_table(out/'regional_observation_equivalence.tsv',checks,list(checks[0]))
    focus={(r['assembly'],r['sequence']) for r in selection if low_retention(r) or ('+' in r['sequence'] and int(r['source_bp'])>=10000000)}
    lines=['# Why chromosome regions were retained or rejected','',
        'This diagnostic observes the current 1 Mb regional assignment rules. It changes no thresholds or assignment conditions. '
        'For each filtering policy and both original/corrected grouping, the observed and unobserved tools must produce '
        'byte-identical PAF and sequence-list files. The production graph is not rebuilt.','',
        '| Source | Filtering | Grouping | Chromosome | Role | Longest interval after overlaps, kb | Candidate runs: pass / length fail / density fail | Final assigned, Mb | Explanation |',
        '|---|---|---|---|---|---:|---|---:|---|']
    for r in rows:
        if (r['assembly'],r['sequence']) not in focus:continue
        lines.append(f"| {r['assembly']} / {r['sequence']} | {r['policy']} | {r['variant']} | {r['chromosome']} | {r['role']} | "
            f"{int(r['post_overlap_max_interval_bp'])/1000:.1f} | {r['locally_accepted_runs']} / {r['length_failed_runs']} / "
            f"{r['density_failed_runs']} | {int(r['final_assigned_span_bp'])/1e6:.3f} | {r['reason']} |")
    lines+=['','## Reading the gates','',
        '- **Overlap resolution:** an interval overlapping a longer interval is removed before regional smoothing.',
        '- **Longest-interval gate:** a secondary run needs at least one premerged interval strictly longer than 1 Mb. '
        'A large total run length alone does not satisfy this existing rule.',
        '- **Density gate:** aligned bases divided by the candidate span must meet the configured minimum. '
        'Span includes the end extension when the candidate owns a scaffold end.',
        '- **Whole-query gate:** local acceptance is provisional; adjusted dominant coverage must exceed its threshold.',
        '- **Small dominant regions:** separate minimum-length filtering can remove pieces of the dominant chromosome.',
        '- Length and density failure counts can overlap. Dominant chromosomes do not undergo the secondary-run gates.', '',
        '[All chromosome summaries](regional_gate_summary.tsv) · [Byte-equivalence checks](regional_observation_equivalence.tsv)', '',
        'Each trace folder contains the exact per-run coordinates, measurements and gate outcomes in '
        '`native_events.tsv.gz` and `fixed_events.tsv.gz`, plus final regions, coordinate audits and build provenance. '
        'Chromosomes with no input mappings cannot be diagnosed as a regional rejection. No presence is inferred from scaffold names.','']
    (out/'regional_gate_report.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='mode',required=True)
    b=sub.add_parser('build');b.add_argument('--source',required=True);b.add_argument('--out',required=True)
    r=sub.add_parser('run')
    for arg in ('prepared','tools','fixed_binary','fixed_info','out'):r.add_argument('--'+arg.replace('_','-'),required=True)
    r.add_argument('--contigs',nargs='+',required=True)
    p=sub.add_parser('report');p.add_argument('--inputs',nargs='+',required=True);p.add_argument('--out',required=True)
    args=vars(parser.parse_args());mode=args.pop('mode');{'build':build,'run':run,'report':report}[mode](**args)


if __name__=='__main__':main()
