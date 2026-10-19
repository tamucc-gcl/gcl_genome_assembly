"""Replay pinned Cactus chromosome assignment, without aligning or changing assemblies."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import shutil
import subprocess
from types import SimpleNamespace
import xml.etree.ElementTree as ET

from pangenome_source_audit import read_table, write_table, query_identity, union_bp

CACTUS_COMMIT = 'd7e5fc1b51327a5a4cdb494c85b0607e3d81a713'
CASES = {'baseline': [], 'whole_relaxed': ['-Q','1'],
         'regional_1Mb': ['-u','1000000','-s'], 'regional_3Mb': ['-u','3000000','-s']}


def identity_aliases(ledger):
    aliases = {}
    for row in ledger:
        name = row['graph_name']
        base, dot, hap = name.rpartition('.')
        for alias in (name, base+'#'+hap if dot and hap.isdigit() else name+'#0'):
            if alias in aliases and aliases[alias] != row['id']:
                raise ValueError('Ambiguous identity alias')
            aliases[alias] = row['id']
    return aliases


def select_queries(rows):
    # Whole queries only: Cactus filtering remains independent of other queries.
    selected = {(r['assembly'],r['sequence']) for r in rows
                if r['status'] != 'represented_in_clip' or '+' in r['sequence']}
    controls = Counter()
    for row in sorted(rows,key=lambda r:(r['assembly'],r['sequence'])):
        if row['chromosome_scale']=='yes' and row['status']=='represented_in_clip' and '+' not in row['sequence']:
            if controls[row['assembly']] < 2:
                selected.add((row['assembly'],row['sequence']))
                controls[row['assembly']] += 1
    return [r for r in rows if (r['assembly'],r['sequence']) in selected]


class LocalJob:
    """File-store adapter for the pinned filter/conversion functions, not a Toil run."""
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True,exist_ok=True)
        self.fileStore = self
        self.memory = 8*1024**3
        self.counter = 0

    def getLocalTempDir(self):
        self.counter += 1
        path = self.root/str(self.counter)
        path.mkdir()
        return str(path)

    def readGlobalFile(self, source, destination):
        Path(destination).symlink_to(Path(source).resolve())
        return destination

    def writeGlobalFile(self, path):
        return str(Path(path).resolve())


def initialize_cactus_tools():
    # Nextflow already supplies the pinned container. Imported Cactus helpers do
    # not run the CLI's setup and otherwise default cactus_call to nested Docker.
    from cactus.shared.common import setupBinaries
    setupBinaries(SimpleNamespace(binariesMode='local',latest=False,workDir=None))


def prepare(sources, identities, paf, gfa, out):
    import cactus
    from cactus.shared.version import cactus_commit
    from cactus.shared.configWrapper import ConfigWrapper
    from cactus.refmap.cactus_graphmap import filter_paf
    from cactus.refmap.cactus_minigraph import minigraph_gfa_from_pansn
    if cactus_commit != CACTUS_COMMIT or importlib.metadata.version('cactus') != '3.2.1':
        raise ValueError('Assignment replay requires the audited Cactus 3.2.1 image/commit')
    initialize_cactus_tools()
    out = Path(out).resolve(); out.mkdir(parents=True,exist_ok=True)
    ledger = read_table(identities)
    references = [r['graph_name'] for r in ledger if r['is_reference'].lower()=='true']
    if len(references) != 1:
        raise ValueError('Expected exactly one audited reference')
    aliases = identity_aliases(ledger)
    selected = select_queries(read_table(sources))
    keys = {(r['assembly'],r['sequence']) for r in selected}
    subset = out/'selected.paf'
    with open(paf,encoding='utf-8') as src, subset.open('w',encoding='utf-8') as dst:
        for line in src:
            query = line.split('\t',1)[0]
            if query_identity(query,aliases) in keys:
                dst.write(line)
    config_path = Path(cactus.__file__).parent/'cactus_progressive_config.xml'
    shutil.copyfile(config_path,out/'effective_config.xml')
    config = ConfigWrapper(ET.parse(config_path).getroot())
    job = LocalJob(out/'adapter_tmp')
    # Call the installed pinned implementation, including its gaffilter overlap pass.
    filtered = filter_paf(job,str(subset),config,reference=references[0])
    shutil.copyfile(filtered,out/'filtered.paf')
    counts, spans = Counter(),defaultdict(list)
    with (out/'filtered.paf').open(encoding='utf-8') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            key = query_identity(fields[0],aliases)
            counts[key] += 1
            spans[key].append((int(fields[2]),int(fields[3])))
    for row in selected:
        key = row['assembly'],row['sequence']
        row['filtered_paf_records'] = counts[key]
        row['filtered_query_covered_bp'] = union_bp(spans[key])
    write_table(out/'selection.tsv',selected,list(selected[0]))
    converted = minigraph_gfa_from_pansn(job,{r['graph_name'] for r in ledger},str(gfa),str(Path(gfa).resolve()))
    with gzip.open(converted,'rb') as src, (out/'cactus.sv.gfa').open('wb') as dst:
        shutil.copyfileobj(src,dst)
    provenance = dict(cactus_package=importlib.metadata.version('cactus'),cactus_commit=cactus_commit,
        binaries_mode='local',
        config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),reference=references[0],
        selected_queries=len(selected),selection='all missing/FULL-only paths, composites, two represented chromosome controls per assembly',
        filter='installed cactus.refmap.cactus_graphmap.filter_paf; whole-query subset',
        graph_conversion='installed cactus.refmap.cactus_minigraph.minigraph_gfa_from_pansn')
    binary=shutil.which('rgfa-split')
    if not binary:raise ValueError('Pinned image lacks rgfa-split')
    provenance['rgfa_split_binary_sha256']=hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    provenance['input_mapping_bytes']=Path(paf).stat().st_size
    provenance['input_graph_bytes']=Path(gfa).stat().st_size
    (out/'preparation.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    shutil.copyfile(identities,out/'identities.tsv')
    # Only the two prepared large intermediates survive in the work directory.
    subset.unlink()
    shutil.rmtree(job.root)


def split_command(gfa,paf,config,reference,contigs,case,prefix,log):
    if case not in CASES:
        raise ValueError('Unknown assignment replay case')
    root = ET.parse(config).getroot()
    split, mapping = root.find('graphmap_split'),root.find('graphmap')
    if split is None or mapping is None:
        raise ValueError('Missing Cactus assignment configuration')
    coverage = split.get('minQueryCoverages','').split()
    thresholds = split.get('minQueryCoverageThresholds','').split()
    if not coverage or len(thresholds)!=len(coverage)-1:
        raise ValueError('Invalid assignment coverage bins')
    cmd = ['rgfa-split','-g',str(gfa),'-p',str(paf),'-b',str(prefix),
           '-Q',split.get('minQueryUniqueness','0'),'-a',split.get('ambiguousName','_AMBIGUOUS_'),
           '-L',str(log),'-r','id='+reference+'|']
    for value in coverage:cmd += ['-n',value]
    for value in thresholds:cmd += ['-T',value]
    if mapping.get('minMAPQ'):cmd += ['-A',mapping.get('minMAPQ')]
    for contig in contigs:
        cmd += ['-c',contig,'-c','id='+reference+'|'+contig]
    # -G is deliberately omitted: it controls graph-file writing, not assignment.
    # Replace Q for the relaxed control rather than passing duplicate flags.
    if case=='whole_relaxed':cmd[cmd.index('-Q')+1] = '1'
    else:cmd += CASES[case]
    return cmd


def fragment(name, aliases, lengths):
    match = re.fullmatch(r'(.*):(\d+)-(\d+)',name)
    base,lo,hi = (match[1],int(match[2])-1,int(match[3])) if match else (name,0,None)
    aid,seq = query_identity(base,aliases)
    if not aid:return None
    key = aid,seq
    if key not in lengths:raise ValueError('Unexpected replay query: '+name)
    hi = lengths[key] if hi is None else hi
    if not 0 <= lo < hi <= lengths[key]:raise ValueError('Invalid replay fragment coordinates: '+name)
    return key,lo,hi


def summarize_split(directory,selection,identities,case,ambiguous):
    selected = read_table(selection)
    lengths = {(r['assembly'],r['sequence']):int(r['source_bp']) for r in selected}
    aliases = identity_aliases(read_table(identities))
    assignments, alignments = defaultdict(list),defaultdict(list)
    regions = []
    for path in sorted(Path(directory).glob('split_*.fa_contigs')):
        chrom = path.name[len('split_'):-len('.fa_contigs')]
        with path.open(encoding='utf-8') as handle:
            for line in handle:
                parsed = fragment(line.strip(),aliases,lengths)
                if parsed is None:continue  # synthetic target nodes
                key,lo,hi = parsed
                assignments[key].append((chrom,lo,hi))
                regions.append(dict(case=case,assembly=key[0],sequence=key[1],chromosome=chrom,
                                    start_0based=lo,end_exclusive=hi))
    for path in sorted(Path(directory).glob('split_*.paf')):
        chrom = path.name[len('split_'):-len('.paf')]
        if chrom==ambiguous:continue
        with path.open(encoding='utf-8') as handle:
            for line in handle:
                fields = line.rstrip().split('\t')
                parsed = fragment(fields[0],aliases,lengths)
                if parsed is None:continue
                key,offset,end = parsed
                lo,hi = offset+int(fields[2]),offset+int(fields[3])
                if not offset <= lo <= hi <= end:raise ValueError('Replay PAF outside assigned fragment')
                alignments[key].append((lo,hi))
    rows = []
    for source in selected:
        key = source['assembly'],source['sequence']
        spans = assignments[key]
        named = [(lo,hi) for chrom,lo,hi in spans if chrom!=ambiguous]
        amb = [(lo,hi) for chrom,lo,hi in spans if chrom==ambiguous]
        rows.append(dict(case=case,assembly=key[0],sequence=key[1],source_bp=source['source_bp'],
            chromosome_scale=source['chromosome_scale'],graph_status=source['status'],
            filtered_query_covered_bp=source['filtered_query_covered_bp'],
            assigned_chromosomes=','.join(sorted({c for c,lo,hi in spans if c!=ambiguous})),
            assigned_span_bp=union_bp(named),ambiguous_span_bp=union_bp(amb),
            surviving_query_covered_bp=union_bp(alignments[key])))
    return rows,regions


def parse_decisions(path,aliases):
    pattern = re.compile(r'^(Query contig is ambiguous: |Assigned (?:ref-)?contig to ([^:]+): )'
        r'(.+?)\s+len=(\d+) cov=([\deE.+-]+) \(vs ([\deE.+-]+)\) '
        r'uf=\s*(infinity|[\deE.+-]+) \(vs ([\deE.+-]+)\)')
    rows=[]
    with open(path,encoding='utf-8') as handle:
        for line in handle:
            match=pattern.match(line)
            if not match:continue
            prefix,chrom,query,length,cov,mincov,uniq,minuniq=match.groups()
            aid,seq=query_identity(query,aliases)
            if not aid:continue
            rejected=prefix.startswith('Query')
            reasons=[]
            if rejected:
                if float(cov)<float(mincov):reasons.append('coverage below threshold')
                if uniq!='infinity' and float(uniq)<float(minuniq):reasons.append('chromosome assignment not unique enough')
                if not reasons:reasons.append('ambiguous; inspect precise log values')
            rows.append(dict(assembly=aid,sequence=seq,decision='ambiguous' if rejected else 'assigned',
                assigned_chromosome=chrom or '',coverage_fraction=cov,min_coverage=mincov,
                uniqueness_ratio=uniq,min_uniqueness=minuniq,reason='; '.join(reasons) if rejected else 'passed whole-scaffold assignment'))
    return rows


def replay(prepared,case,contigs,out):
    prepared,out = Path(prepared),Path(out)
    out.mkdir(parents=True,exist_ok=True)
    temp = out/'split_work';temp.mkdir()
    provenance = json.loads((prepared/'preparation.json').read_text())
    cmd = split_command(prepared/'cactus.sv.gfa',prepared/'filtered.paf',prepared/'effective_config.xml',
        provenance['reference'],contigs,case,temp/'split_',out/'assignment_log.txt')
    (out/'command.json').write_text(json.dumps(cmd,indent=2)+'\n',encoding='utf-8')
    with (out/'tool_stderr.txt').open('w') as err:
        subprocess.run(cmd,check=True,stdout=err,stderr=err)
    ambiguous = ET.parse(prepared/'effective_config.xml').getroot().find('graphmap_split').get('ambiguousName','_AMBIGUOUS_')
    rows,regions = summarize_split(temp,prepared/'selection.tsv',prepared/'identities.tsv',case,ambiguous)
    write_table(out/'assignments.tsv',rows,list(rows[0]))
    write_table(out/'regions.tsv',regions,['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
    decisions=parse_decisions(out/'assignment_log.txt',identity_aliases(read_table(prepared/'identities.tsv')))
    write_table(out/'decisions.tsv',decisions,['assembly','sequence','decision','assigned_chromosome',
        'coverage_fraction','min_coverage','uniqueness_ratio','min_uniqueness','reason'])
    # Keep decisions and coordinates, not duplicate multi-GB split PAFs.
    shutil.rmtree(temp)


def report(directories,out):
    rows = [r for d in directories for r in read_table(Path(d)/'assignments.tsv')]
    if {r['case'] for r in rows} != set(CASES):raise ValueError('Incomplete assignment replay batch')
    by = defaultdict(dict)
    for row in rows:
        key = row['assembly'],row['sequence']
        if row['case'] in by[key]:raise ValueError('Duplicate replay result')
        by[key][row['case']] = row
    if any(set(v)!=set(CASES) for v in by.values()):raise ValueError('Replay query sets differ')
    out = Path(out);out.mkdir(parents=True,exist_ok=True)
    write_table(out/'assignment_comparison.tsv',rows,list(rows[0]))
    lines = ['# Chromosome assignment replay','',
        'These tests use the existing minigraph graph and mappings, the installed pinned Cactus PAF filter, '
        'and rgfa-split. No source assemblies or production graph products are changed.','',
        '| Assembly / source | Case | Chromosomes assigned | Assigned span, bp | Surviving mapped query, bp |',
        '|---|---|---|---:|---:|']
    focus = [v for k,v in sorted(by.items()) if v['baseline']['chromosome_scale']=='yes' and
             v['baseline']['graph_status']!='represented_in_clip']
    for variants in focus:
        for case in CASES:
            r = variants[case]
            lines.append(f"| {r['assembly']} / {r['sequence']} | {case} | {r['assigned_chromosomes'] or 'None'} | {int(r['assigned_span_bp']):,} | {int(r['surviving_query_covered_bp']):,} |")
    baseline = next(Path(d) for d in directories if any(r['case']=='baseline' for r in read_table(Path(d)/'assignments.tsv')))
    decisions={(r['assembly'],r['sequence']):r for r in read_table(baseline/'decisions.tsv')}
    lines += ['', '## Baseline cause recorded by the tool', '',
        '| Assembly / source | Decision | Reason | Coverage / minimum | Uniqueness / minimum |', '|---|---|---|---|---|']
    for variants in focus:
        source=variants['baseline'];decision=decisions.get((source['assembly'],source['sequence']))
        if decision:
            lines.append(f"| {source['assembly']} / {source['sequence']} | {decision['decision']} | {decision['reason']} | {decision['coverage_fraction']} / {decision['min_coverage']} | {decision['uniqueness_ratio']} / {decision['min_uniqueness']} |")
        else:
            explanation='No records survived the pinned PAF filter' if not int(source['filtered_query_covered_bp']) else 'No parsed whole-scaffold decision; inspect assignment log'
            lines.append(f"| {source['assembly']} / {source['sequence']} | Unresolved | {explanation} | — | — |")
    lines += ['', '<details><summary>Other missing source paths by assembly</summary>', '',
        'These are source-query totals for paths missing from the production CLIP graph, '
        'not distinct graph-sequence lengths. All cases use the same selected queries.', '',
        '| Assembly | Case | Filtered query coverage, bp | Surviving query coverage, bp |', '|---|---|---:|---:|']
    totals = defaultdict(lambda:[0,0])
    for r in rows:
        if r['graph_status']!='represented_in_clip':
            total=totals[r['assembly'],r['case']]
            total[0]+=int(r['filtered_query_covered_bp'])
            total[1]+=int(r['surviving_query_covered_bp'])
    for (assembly,case),values in sorted(totals.items()):
        lines.append(f'| {assembly} | {case} | {values[0]:,} | {values[1]:,} |')
    lines += ['', '</details>', '', '## Represented-path controls', '',
        '| Regional case | Represented controls tested | Controls with less surviving mapping than baseline |', '|---|---:|---:|']
    controls=[v for v in by.values() if v['baseline']['graph_status']=='represented_in_clip']
    for case in ('regional_1Mb','regional_3Mb'):
        lower=sum(int(v[case]['surviving_query_covered_bp'])<int(v['baseline']['surviving_query_covered_bp']) for v in controls)
        lines.append(f'| {case} | {len(controls)} | {lower} |')
    lines += ['', 'Any decrease needs inspection in the full comparison table; this count does not assign an acceptable loss threshold.', '',
        '## Decision criteria','',
        'If baseline places the missing source in the ambiguous bin, its assignment_log.txt records '
        'coverage and uniqueness thresholds and the chromosome-specific counts. This establishes '
        'assignment rejection rather than merely inferring it from a missing graph path. '
        'If baseline assigns it normally, investigate the subsequent alignment/graph steps instead.','',
        'The relaxed whole-scaffold control is not a proposed production fix: assigning an entire '
        'scaffold to one chromosome can still discard alignments to the other. Compare surviving '
        'query coverage, not just assigned span. Regional cases must retain both chromosome '
        'portions, keep ambiguous spans explicit and avoid damaging represented chromosome controls. '
        'Assignment success does not establish survival through a later Cactus graph build.','',
        '[Every query/case](assignment_comparison.tsv). Per-case assignment logs and '
        '0-based, half-open region tables are in the sibling case directories. The complete '
        'PAF filter/config provenance is in the preparation directory.','']
    (out/'assignment_report.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode',required=True)
    prep = sub.add_parser('prepare')
    for key in ('sources','identities','paf','gfa','out'):prep.add_argument('--'+key,required=True)
    run = sub.add_parser('replay')
    for key in ('prepared','case','out'):run.add_argument('--'+key,required=True)
    run.add_argument('--contigs',nargs='+',required=True)
    summary = sub.add_parser('report');summary.add_argument('--directories',nargs='+',required=True);summary.add_argument('--out',required=True)
    args = vars(p.parse_args());mode=args.pop('mode')
    if mode=='prepare':prepare(**args)
    elif mode=='replay':replay(**args)
    else:report(**args)


if __name__=='__main__':main()
