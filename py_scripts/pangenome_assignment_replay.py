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

from pangenome_source_audit import read_table, write_table, query_identity, union_bp, low_retention

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
                if r['status'] != 'represented_in_clip' or '+' in r['sequence'] or low_retention(r)}
    controls = Counter()
    for row in sorted(rows,key=lambda r:(r['assembly'],r['sequence'])):
        if (row['chromosome_scale']=='yes' and row['status']=='represented_in_clip'
                and '+' not in row['sequence'] and not low_retention(row)):
            if controls[row['assembly']] < 2:
                selected.add((row['assembly'],row['sequence']))
                controls[row['assembly']] += 1
    return [r for r in rows if (r['assembly'],r['sequence']) in selected]


def merged(intervals):
    result = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if result and lo <= result[-1][1]:
            result[-1] = (result[-1][0], max(hi, result[-1][1]))
        else:
            result.append((lo, hi))
    return result


def interval_difference(left, right):
    """Return A minus B for sorted, disjoint half-open interval lists."""
    result, j = [], 0
    for lo, hi in left:
        while j < len(right) and right[j][1] <= lo:
            j += 1
        pos, k = lo, j
        while k < len(right) and right[k][0] < hi:
            start, end = right[k]
            if start > pos:
                result.append((pos, min(start, hi)))
            pos = max(pos, end)
            if pos >= hi:
                break
            k += 1
        if pos < hi:
            result.append((pos, hi))
    return result


def mapping_intervals(path, aliases, lengths):
    """Read raw PAF/GAF query coordinates once; keep only selected queries."""
    spans, counts = defaultdict(list), Counter()
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) < 12:
                raise ValueError('Malformed coverage mapping record')
            key = query_identity(fields[0], aliases)
            if key not in lengths:
                continue
            length, lo, hi = map(int, fields[1:4])
            if length != lengths[key] or not 0 <= lo <= hi <= length:
                raise ValueError('Coverage mapping coordinates disagree with source: '+str(key))
            counts[key] += 1
            spans[key].append((lo, hi))
    return {key: merged(spans[key]) for key in lengths}, counts


def write_intervals(path, stages):
    write_table(path, (dict(stage=stage, assembly=key[0], sequence=key[1],
                           start_0based=lo, end_exclusive=hi)
                       for stage, by_key in stages.items() for key, spans in sorted(by_key.items())
                       for lo, hi in spans),
                ['stage','assembly','sequence','start_0based','end_exclusive'])


def bin_coverage(spans, length, size=1000000):
    values = [0] * ((length+size-1)//size)
    for lo, hi in spans:
        while lo < hi:
            idx = lo//size
            end = min(hi, (idx+1)*size)
            values[idx] += end-lo
            lo = end
    return values


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


def prepare(sources, identities, paf, gfa, gaf, out):
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
    lengths = {(r['assembly'],r['sequence']):int(r['source_bp']) for r in selected}
    raw_gaf, _ = mapping_intervals(gaf, aliases, lengths)
    raw_paf, _ = mapping_intervals(paf, aliases, lengths)
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
    spans, counts = mapping_intervals(out/'filtered.paf', aliases, lengths)
    write_intervals(out/'coverage_intervals.tsv.gz',
                    dict(raw_gaf=raw_gaf, raw_paf=raw_paf, filtered_paf=spans))
    for row in selected:
        key = row['assembly'],row['sequence']
        row['filtered_paf_records'] = counts[key]
        row['filtered_query_covered_bp'] = union_bp(spans[key])
        row['raw_gaf_query_covered_bp'] = union_bp(raw_gaf[key])
        row['raw_paf_query_covered_bp'] = union_bp(raw_paf[key])
    write_table(out/'selection.tsv',selected,list(selected[0]))
    converted = minigraph_gfa_from_pansn(job,{r['graph_name'] for r in ledger},str(gfa),str(Path(gfa).resolve()))
    with gzip.open(converted,'rb') as src, (out/'cactus.sv.gfa').open('wb') as dst:
        shutil.copyfileobj(src,dst)
    provenance = dict(cactus_package=importlib.metadata.version('cactus'),cactus_commit=cactus_commit,
        binaries_mode='local',
        config_sha256=hashlib.sha256(config_path.read_bytes()).hexdigest(),reference=references[0],
        selected_queries=len(selected),selection='all missing/FULL-only paths, composites, chromosome-scale paths below 80% retained length, two represented chromosome controls per assembly',
        coverage_coordinates='0-based half-open query intervals; union coverage, not alignment sum',
        input_gaf_bytes=Path(gaf).stat().st_size,
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


def summarize_split(directory,selection,identities,case,ambiguous,interval_output=None):
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
    if interval_output is not None:
        write_intervals(interval_output, {case:{key:merged(alignments[key]) for key in lengths}})
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
        for field in ('raw_gaf_query_covered_bp','raw_paf_query_covered_bp',
                      'clip_retained_path_bp','retention_review'):
            rows[-1][field] = source.get(field, '')
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
    rows,regions = summarize_split(temp,prepared/'selection.tsv',prepared/'identities.tsv',case,ambiguous,
                                   out/'surviving_intervals.tsv.gz')
    write_table(out/'assignments.tsv',rows,list(rows[0]))
    write_table(out/'regions.tsv',regions,['case','assembly','sequence','chromosome','start_0based','end_exclusive'])
    decisions=parse_decisions(out/'assignment_log.txt',identity_aliases(read_table(prepared/'identities.tsv')))
    write_table(out/'decisions.tsv',decisions,['assembly','sequence','decision','assigned_chromosome',
        'coverage_fraction','min_coverage','uniqueness_ratio','min_uniqueness','reason'])
    # Keep decisions and coordinates, not duplicate multi-GB split PAFs.
    shutil.rmtree(temp)


def coverage_report(directories, coverage, selection, out):
    sources = read_table(selection)
    lengths = {(r['assembly'],r['sequence']):int(r['source_bp']) for r in sources}
    stages = {name:defaultdict(list) for name in ('raw_gaf','raw_paf','filtered_paf',*CASES)}
    regions = defaultdict(list)
    paths = [Path(coverage)] + [Path(d)/'surviving_intervals.tsv.gz' for d in directories]
    for path in paths:
        for row in read_table(path):
            key = row['assembly'], row['sequence']
            stage = row['stage']
            lo, hi = int(row['start_0based']), int(row['end_exclusive'])
            if key not in lengths or stage not in stages or not 0 <= lo < hi <= lengths[key]:
                raise ValueError('Invalid coverage trace interval')
            stages[stage][key].append((lo,hi))
    for d in directories:
        for r in read_table(Path(d)/'regions.tsv'):
            regions[r['case'],r['assembly'],r['sequence']].append(r)
    for values in stages.values():
        for key in values:
            values[key] = merged(values[key])
    summary, bins, changes = [], [], []
    previous = dict(raw_gaf=None, raw_paf='raw_gaf', filtered_paf='raw_paf',
                    **{case:'filtered_paf' for case in CASES})
    for source in sources:
        key = source['assembly'],source['sequence']
        for stage, values in stages.items():
            spans = values[key]
            prior = previous[stage]
            lost = interval_difference(stages[prior][key],spans) if prior else []
            gained = interval_difference(spans,stages[prior][key]) if prior else []
            summary.append(dict(assembly=key[0],sequence=key[1],stage=stage,source_bp=lengths[key],
                covered_bp=union_bp(spans),covered_fraction=union_bp(spans)/lengths[key],
                compared_with=prior or '',lost_covered_bp=union_bp(lost) if prior else '',
                gained_covered_bp=union_bp(gained) if prior else ''))
            for direction, intervals in (('lost',lost),('gained',gained)):
                changes.extend(dict(assembly=key[0],sequence=key[1],stage=stage,compared_with=prior,
                    change=direction,start_0based=lo,end_exclusive=hi) for lo,hi in intervals)
            for i, covered in enumerate(bin_coverage(spans,lengths[key])):
                lo, hi = i*1000000, min((i+1)*1000000,lengths[key])
                chroms = sorted({r['chromosome'] for r in regions[(stage,*key)]
                    if int(r['start_0based']) < hi and int(r['end_exclusive']) > lo})
                bins.append(dict(assembly=key[0],sequence=key[1],stage=stage,start_0based=lo,
                    end_exclusive=hi,covered_bp=covered,covered_fraction=covered/(hi-lo),
                    assigned_chromosomes=','.join(chroms)))
    write_table(out/'coverage_trace.tsv',summary,list(summary[0]))
    write_table(out/'coverage_bins.tsv',bins,list(bins[0]))
    write_table(out/'coverage_changes.tsv.gz',changes,
        ['assembly','sequence','stage','compared_with','change','start_0based','end_exclusive'])
    # Preserve the merged source-coordinate intervals for reproduction without huge mappings.
    write_intervals(out/'coverage_intervals.tsv.gz',stages)
    by = {(r['assembly'],r['sequence'],r['stage']):r for r in summary}
    lines = ['', '## Where source coverage is lost', '',
        'Union coverage of query coordinates, in Mb. Zero means the supplied mapping file '
        'was read and contained no coverage for this selected source. Reference query mapping '
        'can be absent because reference paths are extracted from graph tags; zero is not evidence of reference loss.', '',
        '| Assembly / source | Input | Raw GAF | Raw PAF | Filtered PAF | Baseline assigned mapping | CLIP path length |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in sources:
        key = r['assembly'],r['sequence']
        if '+' not in key[1] and not (r['chromosome_scale']=='yes' and
                (r['status']!='represented_in_clip' or low_retention(r))):
            continue
        nums = [int(r['source_bp'])] + [by[(*key,s)]['covered_bp'] for s in
            ('raw_gaf','raw_paf','filtered_paf','baseline')]
        values = ' | '.join(f'{n/1e6:.3f}' for n in nums)
        clip = r.get('clip_retained_path_bp','')
        clip_text = f'{int(clip)/1e6:.3f}' if clip != '' else 'Unavailable'
        lines.append(f'| {key[0]} / {key[1]} | {values} | {clip_text} |')
    lines += ['', 'Loss and gain are coordinate set differences between consecutive mapping stages; '
        'each assignment case is compared directly with the filtered PAF. GAF-to-PAF conversion '
        'and filtering can alter coverage, so net totals alone can hide relocation. '
        'No coverage threshold changes the production assignment policy.', '',
        'CLIP path length is a traversal-length measurement, not query-coordinate coverage. '
        'It must not be subtracted from mapping coverage to claim precise clipping losses. '
        'FULL W-coordinate spans likewise do not certify complete source coverage.', '',
        '[All source/stage totals and losses/gains](coverage_trace.tsv) · '
        '[1 Mb source-coordinate coverage and assignment context](coverage_bins.tsv) · '
        '[Exact lost/gained query intervals](coverage_changes.tsv.gz) · '
        '[Merged covered intervals](coverage_intervals.tsv.gz)', '',
        'Coordinates are 0-based, half-open. Bin chromosome labels describe overlapping assigned '
        'regions, including ambiguous regions; they do not imply every mapped base has that label. '
        'The merged intervals and changes resolve boundaries more precisely than the 1 Mb summaries.', '']
    return lines


def report(directories,out,coverage=None,selection=None):
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
    focus = [v for k,v in sorted(by.items()) if '+' in k[1] or
             (v['baseline']['chromosome_scale']=='yes' and
              (v['baseline']['graph_status']!='represented_in_clip' or
               v['baseline'].get('retention_review')=='yes'))]
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
    controls=[v for key,v in by.items() if v['baseline']['graph_status']=='represented_in_clip'
              and '+' not in key[1] and v['baseline'].get('retention_review')!='yes']
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
    if coverage is not None and selection is not None:
        lines += coverage_report(directories, coverage, selection, out)
    (out/'assignment_report.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode',required=True)
    prep = sub.add_parser('prepare')
    for key in ('sources','identities','paf','gfa','gaf','out'):prep.add_argument('--'+key,required=True)
    run = sub.add_parser('replay')
    for key in ('prepared','case','out'):run.add_argument('--'+key,required=True)
    run.add_argument('--contigs',nargs='+',required=True)
    summary = sub.add_parser('report');summary.add_argument('--directories',nargs='+',required=True);summary.add_argument('--out',required=True)
    summary.add_argument('--coverage',required=True);summary.add_argument('--selection',required=True)
    args = vars(p.parse_args());mode=args.pop('mode')
    if mode=='prepare':prepare(**args)
    elif mode=='replay':replay(**args)
    else:report(**args)


if __name__=='__main__':main()
