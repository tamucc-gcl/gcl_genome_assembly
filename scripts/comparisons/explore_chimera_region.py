#!/usr/bin/env python3
"""Read-only regional diagnostics from existing Nextflow work files. Run remotely."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
from collections import Counter
from pathlib import Path


def table(path):
    with open(path) as f:
        return list(csv.DictReader((s for s in f if s.strip() and not s.startswith('#')), delimiter='\t'))


def write_table(path, fields, rows):
    with open(path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def one(paths, label):
    paths = list(paths)
    if len(paths) != 1:
        raise ValueError(f'Expected one {label}, found {paths}')
    return paths[0]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def assessment_fasta(work_dirs, expected, explicit=None):
    """Resolve staged FASTA by content, not Nextflow's staging filename.

    Search only task roots and one directory level, never the whole work tree.
    Identical copies/links are interchangeable; a wrong-stage FASTA is rejected.
    """
    if explicit is not None:
        if sha(explicit) != expected:
            raise ValueError('Explicit assessment FASTA checksum mismatch')
        return explicit
    candidates, inventory = [], []
    for root in work_dirs:
        for path in sorted(root.iterdir()):
            inventory.append(str(path))
            if path.is_file():
                candidates.append(path)
            elif path.is_dir() and not path.is_symlink():
                candidates.extend(p for p in sorted(path.iterdir()) if p.is_file())
    seen = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        with open(path, 'rb') as f:
            is_fasta = f.read(1) == b'>'
        if is_fasta and sha(path) == expected:
            return path
    raise ValueError('No staged FASTA matches the recorded assessment SHA256. '
                     'Retain the task inputs or supply --assessment-fasta with the exact '
                     'pre-finishing FASTA. Task entries: '+', '.join(inventory))


def local_contacts(counts, center, flank, distance_lo, distance_hi):
    """Equal-distance comparison; return counts AND available bin-pair counts.

    Unbalanced exploratory statistic, not a test or calibrated cut threshold.
    Every bin pair is counted once; self-bin contacts are excluded.
    """
    totals = Counter()
    opportunities = Counter()
    for x in range(center-flank, center+flank):
        for y in range(x+1, center+flank):
            if not distance_lo <= y-x < distance_hi:
                continue
            category = 'cross' if x < center <= y else 'left' if y < center else 'right'
            totals[category] += counts.get((x, y), 0)
            opportunities[category] += 1
    return totals, opportunities


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--joins-work', type=Path, required=True)
    p.add_argument('--evidence-work', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--assessment-fasta', type=Path,
                   help='optional exact pre-finishing FASTA; SHA256 must match the call table')
    p.add_argument('--comparison-assembly', action='append', default=[])
    p.add_argument('--log', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--scaffold', required=True)
    p.add_argument('--start', type=int, required=True, help='0-based region start')
    p.add_argument('--end', type=int, required=True, help='exclusive region end')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not 0 <= a.start < a.end:
        p.error('Require 0 <= start < end')
    a.out.mkdir(parents=True, exist_ok=False)
    calls = a.joins_work / (a.assembly+'.chimeric_joins.tsv')
    selected = [r for r in table(calls) if r['scaffold'] == a.scaffold]
    if not selected:
        raise ValueError('No matching candidate in this work directory')
    digests = {r['assembly_sha256'] for r in selected}
    if len(digests) != 1 or not re.fullmatch('[0-9a-f]{64}', next(iter(digests))):
        raise ValueError('Candidate rows do not identify one valid assessment SHA256')
    digest = next(iter(digests))
    fasta = assessment_fasta([a.joins_work, a.evidence_work], digest, a.assessment_fasta)
    evidence_calls = one(list(a.evidence_work.glob('*.review_joins.tsv'))+
                         list(a.evidence_work.glob('*.chimeric_joins.tsv')), 'evidence call table')
    evidence_rows = [r for r in table(evidence_calls) if r['scaffold'] == a.scaffold]
    if not evidence_rows or any(r['assembly_sha256'] != digest for r in evidence_rows):
        raise ValueError('Evidence work directory refers to a different assessment FASTA')
    seq, active = [], False
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                active = line[1:].split()[0] == a.scaffold
            elif active:
                seq.append(line.strip())
    seq = ''.join(seq).upper()
    if a.end > len(seq):
        raise ValueError('Requested region exceeds scaffold length')
    shutil.copy2(calls, a.out/'native_calls.tsv')
    gaps = [dict(start=m.start(), end=m.end(), length=m.end()-m.start())
            for m in re.finditer('N+', seq) if m.start() < a.end and m.end() > a.start]
    write_table(a.out/'sequence_gaps.tsv', ['start', 'end', 'length'], gaps)
    agp = one(a.joins_work.glob('*.agp'), 'last-round AGP')
    with open(agp) as source, open(a.out/'region.agp', 'w') as dest:
        for line in source:
            f = line.split()
            if f and not line.startswith('#') and f[0] == a.scaffold and int(f[1])-1 < a.end and int(f[2]) > a.start:
                dest.write(line)
    paf = one(a.joins_work.glob('*.ref.paf*'), 'reference PAF')
    opener = gzip.open if paf.suffix == '.gz' else open
    targets = set()
    fields = ['query', 'qlen', 'qstart', 'qend', 'strand', 'target', 'tlen', 'tstart', 'tend', 'matches', 'block', 'mapq', 'tags']
    regional = []
    with opener(paf, 'rt') as f:
        for line in f:
            row = line.rstrip().split('\t')
            if len(row) < 12:
                raise ValueError('Malformed PAF')
            if row[0] == a.scaffold and int(row[2]) < a.end and int(row[3]) > a.start:
                targets.add(row[5])
                regional.append(dict(zip(fields, row[:12]+[';'.join(row[12:])])))
    write_table(a.out/'regional_alignments.tsv', fields, regional)
    # Keep all mappings to these reference targets: complementary/duplicated arms
    # must be judged using intervals, not merely assigned chromosome names.
    with opener(paf, 'rt') as source, open(a.out/'target_context.paf', 'w') as dest:
        for line in source:
            if line.split('\t')[5] in targets:
                dest.write(line)
    for path in a.joins_work.glob('*name_map.tsv'):
        shutil.copy2(path, a.out/path.name)
    for assembly in a.comparison_assembly:
        path = a.results/'assembly/harmonization'/(assembly+'.ref.paf.gz')
        with gzip.open(path, 'rt') as source, open(a.out/(assembly+'.target_context.paf'), 'w') as dest:
            for line in source:
                if line.split('\t')[5] in targets:
                    dest.write(line)

    # Existing pairs are 1-based; convert before binning. Read the full scaffold
    # once, retain only the requested region plus full two-sided flank windows.
    pairs = a.evidence_work/(a.assembly+'.'+a.scaffold+'.pairs')
    resolution, flank = 100000, 20
    first = max(flank, (a.start+resolution-1)//resolution)
    last = min(len(seq)//resolution-flank, a.end//resolution)
    counts, marginals = Counter(), Counter()
    if not pairs.is_file():
        raise ValueError('Expected existing projected pairs are missing: '+str(pairs))
    with open(pairs) as f:
        for line in f:
            if line.startswith('#'):
                continue
            r = line.split()
            if r[1] != a.scaffold or r[3] != a.scaffold:
                raise ValueError('Unexpected scaffold in projected pairs')
            x, y = sorted(((int(r[2])-1)//resolution, (int(r[4])-1)//resolution))
            marginals[x] += 1
            marginals[y] += 1
            if first-flank <= x <= y < last+flank:
                counts[x, y] += 1
    profiles = []
    for center in range(first, last+1):
        for low, high in ((1, 5), (5, 10), (10, 20)):
            total, possible = local_contacts(counts, center, flank, low, high)
            row = dict(position=center*resolution, distance_min=low*resolution, distance_max=high*resolution)
            for label in ('cross', 'left', 'right'):
                row[label+'_contacts'] = total[label]
                row[label+'_bin_pairs'] = possible[label]
                row[label+'_mean'] = total[label]/possible[label] if possible[label] else ''
            within_n = possible['left']+possible['right']
            within = (total['left']+total['right'])/within_n if within_n else 0
            row['cross_over_within'] = row['cross_mean']/within if within and row['cross_mean'] != '' else ''
            profiles.append(row)
    if not profiles:
        raise ValueError('Region has no full-window Hi-C positions')
    write_table(a.out/'local_hic.tsv', list(profiles[0]), profiles)
    write_table(a.out/'hic_bin_coverage.tsv', ['start', 'end', 'read_ends'],
                [dict(start=b*resolution, end=min((b+1)*resolution,len(seq)), read_ends=marginals[b])
                 for b in range(first-flank, last+flank)])
    write_table(a.out/'hic_matrix.tsv', ['bin1_start', 'bin2_start', 'contacts'],
                [dict(bin1_start=x*resolution, bin2_start=y*resolution, contacts=n)
                 for (x,y),n in sorted(counts.items())])
    # Reports belong to earlier coordinates: copy as labelled context, never
    # interpret their scaffold names as current positions or run their commands.
    inventory = []
    for stage in ('contig', 'scaffold'):
        parent = a.results/'assembly'/stage/'misassembly_correction'
        for path in sorted(parent.glob(a.assembly+'*')):
            if not path.is_file():
                continue
            inventory.append(dict(stage=stage, path=str(path.resolve()), bytes=path.stat().st_size))
            if path.suffix in ('.bed', '.txt'):
                dest = a.out/'inspector_context'/stage
                dest.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest/path.name)
    write_table(a.out/'inspector_inventory.tsv', ['stage', 'path', 'bytes'], inventory)
    # Resolve the cached Inspector tasks from this invocation; collect their exact
    # input/output links and command text, without running any task shell script.
    work_links = []
    with open(a.log) as f:
        for line in f:
            if 'CORRECT_MISASSEMBLIES' not in line or '('+a.assembly+')' not in line:
                continue
            match = re.search(r'\[([0-9a-f]{2})/([0-9a-f]{6})\].*(?:Cached|Submitted) process', line)
            if not match:
                continue
            task = one((a.work_root/match[1]).glob(match[2]+'*'), 'Inspector task directory')
            dest = a.out/'inspector_context'/('task_'+match[1]+'_'+match[2])
            dest.mkdir(parents=True, exist_ok=True)
            shutil.copy2(task/'.command.sh', dest/'command.txt')
            for path in sorted(task.iterdir()):
                if path.is_file() and path.suffix in ('.fa', '.fasta', '.fq', '.fastq', '.gz', '.bam'):
                    work_links.append(dict(task=str(task), name=path.name,
                                           resolved=str(path.resolve()), bytes=path.stat().st_size))
    write_table(a.out/'inspector_work_inputs.tsv', ['task', 'name', 'resolved', 'bytes'], work_links)
    older = a.results/'assembly/chimeras/older_joins'/(a.assembly+'.older_join_audit.tsv')
    if older.is_file():
        records = [r for r in table(older) if r['current_object'] == a.scaffold
                   and r['current_lo'] and int(r['current_lo']) < a.end and int(r['current_hi']) > a.start]
        if records:
            write_table(a.out/'regional_older_joins.tsv', list(records[0]), records)
    for path in a.evidence_work.glob(a.assembly+'.'+a.scaffold+'*.chimera_evidence.*'):
        shutil.copy2(path, a.out/path.name)
    (a.out/'provenance.json').write_text(json.dumps(dict(
        assembly=a.assembly, scaffold=a.scaffold, region=[a.start,a.end],
        fasta=str(fasta.resolve()), sha256=digest, scaffold_length=len(seq),
        paf=str(paf.resolve()), pairs=str(pairs.resolve()), targets=sorted(targets),
        joins_work=str(a.joins_work.resolve()), evidence_work=str(a.evidence_work.resolve()),
        coordinate_system='0-based half-open except original AGP and PAF/SAM tool conventions',
        resolution=resolution, flank_bp=flank*resolution,
        note='Exploratory unbalanced Hi-C ratios, not a cut test. HiFi support not yet adjudicated.'), indent=2)+'\n')
    print('Diagnostics written to', a.out)


if __name__ == '__main__':
    main()
