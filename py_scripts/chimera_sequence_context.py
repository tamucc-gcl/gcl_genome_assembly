#!/usr/bin/env python3
"""Published-tool diagnostic orchestration. Never emits permission to cut."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import re


def eligible(rows):
    return [r for r in rows if r.get('chromosome_member') == 'yes'
            and r.get('candidate_verdict') in ('REVIEW', 'BREAK_CANDIDATE')
            and (r.get('callable') == 'yes' or r.get('evidence_only') == 'yes')]


def spans_interval(fields, lo, hi, flank=1000):
    """A single gapped PAF record brackets an interval; NOT a proof of continuity."""
    return int(fields[2]) <= lo - flank and int(fields[3]) >= hi + flank


def narrow_spanning_molecules(path, lo, hi, flank=1000):
    """Primary MAPQ30 molecules with uninterrupted aligned sequence over both anchors.

    Reference deletions/skips over 50 bp veto a span. Broad bracketing alone is insufficient.
    This is positive local support only; zero is not informative absence.
    """
    molecules = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.split('\t')
            if len(f) < 11 or int(f[1]) & (4 | 256 | 2048) or int(f[4]) < 30:
                continue
            start = int(f[3])-1
            position, blocks = start, []
            bad = False
            for length, operation in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
                length = int(length)
                if operation in 'M=X':
                    blocks.append((position, position+length)); position += length
                elif operation in 'DN':
                    if length > 50 and position < hi+flank and position+length > lo-flank:
                        bad = True
                    position += length
            left = sum(max(0, min(end, lo)-max(begin, lo-flank)) for begin, end in blocks)
            right = sum(max(0, min(end, hi+flank)-max(begin, hi)) for begin, end in blocks)
            if not bad and start <= lo-flank and position >= hi+flank and left >= flank and right >= flank:
                molecules.add(f[0])
    return len(molecules)


def run(args, output=None):
    if output:
        with open(output, 'w') as handle:
            subprocess.run([str(x) for x in args], stdout=handle, check=True)
    else:
        subprocess.run([str(x) for x in args], check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--calls', required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--peers', required=True, help='JSON list of id/sample/path records')
    p.add_argument('--reads', default='')
    p.add_argument('--motif', required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir()
    with open(a.calls) as handle:
        rows = eligible(list(csv.DictReader((x for x in handle if not x.startswith('#')), delimiter='\t')))
    if not rows:
        (out/'report.md').write_text('# Chimera sequence context\n\nNo eligible chromosome-scale candidate intervals.\n')
        return
    digest = hashlib.sha256()
    with open(a.fasta, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b''):
            digest.update(chunk)
    if any(r.get('assembly_sha256') != digest.hexdigest() for r in rows):
        raise ValueError('Candidate assessment checksum does not match FASTA')
    # Index a local symlink; never write an index next to a published input.
    local = out/'assessment.fa'
    local.symlink_to(Path(a.fasta).resolve())
    run(['samtools', 'faidx', local])
    sizes = {f[0]: int(f[1]) for f in (x.split('\t') for x in Path(str(local)+'.fai').read_text().splitlines())}
    if sum(sizes.values()) >= 8_000_000_000:
        raise ValueError('Assessment exceeds diagnostic single-index budget')
    queries = out/'intervals.fa'
    intervals = {}
    decision_measurements = {}
    with queries.open('w') as combined:
        for i, row in enumerate(rows):
            lo, hi = int(row['transition_lo']), int(row['transition_hi'])
            chrom = row['scaffold']
            if not 0 <= lo <= hi <= sizes[chrom]:
                raise ValueError('Invalid assessment interval')
            start, end = max(0, lo-100000), min(sizes[chrom], hi+100000)
            key = 'candidate_%d' % (i+1)
            fa = out/(key+'.fa')
            run(['samtools', 'faidx', local, '%s:%d-%d' % (chrom, start+1, end)], fa)
            seq = ''.join(fa.read_text().splitlines()[1:])
            fa.write_text('>'+key+'\n'+seq+'\n')
            combined.write(fa.read_text())
            intervals[key] = dict(scaffold=chrom, start=start, end=end, lo=lo, hi=hi)
            decision_id = hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()[:20]
            intervals[key]['decision_id'] = decision_id
            decision_measurements[decision_id] = dict(assessment_sha256=digest.hexdigest(),
                                                      coordinate_stage='pre_finishing')
    run(['tidk', 'search', '--string', a.motif, '--window', '100', '--output', 'intervals', '--dir', out/'tidk', queries])
    peers = json.loads(a.peers)
    summaries = {key: dict(same=set(), other=set()) for key in intervals}
    with (out/'peer_context.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['candidate', 'peer', 'sample', 'relationship', 'single_record_brackets_interval', 'mapq', 'target', 'target_start', 'target_end', 'target_length'])
        for n, peer in enumerate(peers):
            peer_local = out/('peer_%d.fa' % n)
            peer_local.symlink_to(Path(peer['path']).resolve())
            run(['samtools', 'faidx', peer_local])
            length = sum(int(line.split('\t')[1]) for line in Path(str(peer_local)+'.fai').read_text().splitlines())
            if length >= 8_000_000_000:
                raise ValueError('Peer exceeds diagnostic single-index budget')
            paf = out/('peer_%d.paf' % n)
            run(['minimap2', '-x', 'asm5', '-I', '8G', '-c', '--secondary=yes', '-N', '20', '-t', a.threads, peer['path'], queries], paf)
            # Report every alignment, including ambiguity; a missing row is not negative support.
            for line in paf.read_text().splitlines():
                f = line.split('\t')
                interval = intervals[f[0]]
                bracket = spans_interval(f, interval['lo']-interval['start'], interval['hi']-interval['start'])
                if bracket and int(f[11]) >= 20:
                    relationship = 'same' if peer['sample'] == a.sample else 'other'
                    summaries[f[0]][relationship].add(peer['sample'])
                writer.writerow([f[0], peer['id'], peer['sample'], 'same_individual' if peer['sample'] == a.sample else 'other_individual', 'yes' if bracket else 'no', f[11], f[5], f[7], f[8], f[6]])
            peer_local.unlink()
    if a.reads:
        sam = out/'mapping.sam'
        run(['minimap2', '-ax', 'map-hifi', '-I', '8G', '--secondary=yes', '-N', '20', '-t', a.threads, local, a.reads], sam)
        bam = out/'hifi.bam'
        run(['samtools', 'sort', '-@', '2', '-m', '2G', '-o', bam, sam])
        sam.unlink()
        run(['samtools', 'index', bam])
        run(['samtools', 'quickcheck', bam])
        for key, interval in intervals.items():
            region = '%s:%d-%d' % (interval['scaffold'], interval['start']+1, interval['end'])
            run(['samtools', 'view', '-h', bam, region], out/(key+'.sam'))
            decision_measurements[interval['decision_id']]['hifi_spanning_molecules'] = narrow_spanning_molecules(
                out/(key+'.sam'), interval['lo'], interval['hi'])
            run(['samtools', 'depth', '-aa', '-q', '0', '-Q', '20', '-G', '0xF04', '-r', region, bam], out/(key+'.depth.tsv'))
    (out/'decision_measurements.json').write_text(json.dumps(decision_measurements, indent=2)+'\n')
    for tool in ('samtools', 'minimap2', 'tidk'):
        run([tool, '--version'], out/(tool+'.version.txt'))
    local.unlink()
    (out/'provenance.json').write_text(json.dumps(dict(assembly=a.assembly, sample=a.sample, sha256=digest.hexdigest(), coordinate_stage='pre_finishing', intervals=intervals, peers=peers, hifi=bool(a.reads)), indent=2)+'\n')
    (out/'report.md').write_text(
        '# Chimera sequence context: '+a.assembly+'\n\n'
        +str(len(intervals))+' candidate intervals; '+str(len(peers))+' same-species peer assemblies.\n\n'
        'HiFi mapping: '+('available' if a.reads else 'not requested or no HiFi reads')+'.\n\n'
        'See peer_context.tsv for direct interval alignments and provenance.json for coordinate offsets. '
        'A single gapped alignment bracketing an interval is contextual support, not proof of basewise continuity. '
        'Same-individual haplotypes are not independent votes. No alignment is not evidence of a technical error.\n\n'
        'tidk windows describe motif counts, not repeat-tract lengths or chromosome fusions. '
        'HiFi evidence requires inspection of anchors, repeats and alternative placements. '
        'No diagnostic here authorizes a cut; candidate and cutting decisions remain separate.\n')
    with (out/'report.md').open('a') as report:
        report.write('\n<details>\n<summary>Interval-level peer context</summary>\n\n'
                     'Counts below require MAPQ >=20 and a single alignment with 1-kb flanks. '
                     'They count distinct individuals, not alignment records, and are not cutting votes.\n\n'
                     '| Candidate | Scaffold interval (0-based, half-open) | Same individual observed | Other individuals observed |\n'
                     '|---|---|---|---:|\n')
        for key, interval in intervals.items():
            summary = summaries[key]
            report.write('| %s | %s:%d-%d | %s | %d |\n' %
                         (key, interval['scaffold'], interval['lo'], interval['hi'],
                          'yes' if summary['same'] else 'not observed', len(summary['other'])))
        report.write('\nZero or not observed means no qualifying record, not demonstrated discontinuity.\n\n</details>\n')


if __name__ == '__main__':
    main()
