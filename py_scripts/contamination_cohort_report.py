#!/usr/bin/env python3
"""Present FCS screening and optional exhaustive BlobTools evidence from this run."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from urllib.parse import quote


def rows(path):
    with Path(path).open(encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, records, fields):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def fcs_rows(path):
    lines = Path(path).read_text().splitlines()
    header = next((x for x in lines if x.startswith('#seq_id\t') or x.startswith('#seq-id\t')), None)
    if header is None:
        raise ValueError(f'FCS action header missing: {path}')
    return list(csv.DictReader([header.lstrip('#')] + [x for x in lines if x and not x.startswith('#')], delimiter='\t'))


def adaptor_rows(path):
    lines = Path(path).read_text().splitlines()
    return [dict(zip(['accession','length','action','range','name'],line.split('\t',4)))
            for line in lines if line and not line.startswith('#') and
            line.split('\t')[0] != 'accession']


def unique(records, label):
    result = {}
    for row in records:
        name = str(row['assembly'])
        if name in result:
            raise ValueError(f'Duplicate {label} assembly: {name}')
        result[name] = row
    return result


def build(spec, out):
    out.mkdir(parents=True, exist_ok=True)
    detail = out / 'report'
    detail.mkdir(exist_ok=True)
    expected = sorted(spec['assemblies'])
    if len(expected) != len(set(expected)):
        raise ValueError('Duplicate expected assemblies')
    fcs = unique(spec.get('fcs', []), 'FCS')
    blob = unique(spec.get('blobtools', []), 'BlobTools')
    metrics = unique(spec.get('fcs_metrics', []), 'FCS metrics')
    adapters = unique(spec.get('adaptor_reports', []), 'FCS-adaptor')
    for label, enabled, records in [('FCS', spec['fcs_enabled'], fcs), ('BlobTools', spec['blobtools_enabled'], blob)]:
        if enabled and set(records) != set(expected):
            raise ValueError(f'{label} current-run membership differs: expected {expected}, received {sorted(records)}')
        if not enabled and records:
            raise ValueError(f'{label} supplied despite being disabled')
    for label, records in [('FCS metrics',metrics),('FCS-adaptor',adapters)]:
        if records and (not spec['fcs_enabled'] or set(records) != set(expected)):
            raise ValueError(label+' membership differs from the FCS cohort')
    summary, actions, adapter_actions, intervals, sequences, fcs_sequences = [], [], [], [], [], []
    final = spec['stage'] == 'final'
    lines = ['# Contamination evidence across assemblies', '',
        f"Assembly stage: **{spec['stage']}**. FCS: **{'run' if spec['fcs_enabled'] else 'not run'}**. "
        f"Optional BlobTools searches: **{'run' if spec['blobtools_enabled'] else 'not run'}**.", '']
    if final:
        lines += ['**Final assembly audit:** adaptor and GX screening use the exact finalized FASTA after gap filling, '
                  'telomere extension and final naming/reorientation. No sequence is edited by this audit. '
                  'Any new cleaning recommendation must be resolved through tracked cleanup and re-screening; '
                  'a successful task exit does not mean the assembly has no findings.', '']
    else:
        lines += ['FCS-adaptor screens the stage input and its cleaning step removes terminal adaptors, splits at '
                  'internal adaptor spans, and excludes adaptor-dominated sequences. GX then screens that cleaned input '
                  'and performs its own cleaning. FCS plot coordinates refer to the exact GX input, **before GX cleanup**, '
                  'so sequences subsequently removed by GX are visible.', '']
    lines += ['FCS action coordinates are 1-based inclusive. Optional BlobTools evidence uses the post-FCS assembly '
              'and 0-based half-open coordinates. Do not transfer intervals between these sequence versions. '
              'FCS action counts are recommendations from the reports, not a measurement of bases removed.', '',
              '| Assembly | Adaptor actions | GX actions | GX review actions | Audit status | FCS results |',
              '|---|---:|---:|---:|---|---|']
    for name in expected:
        row = dict(assembly=name,stage=spec['stage'],fcs_status='not run',adaptor_status='not run',
            adaptor_actions='',fcs_actions='',fcs_review_actions='',audit_status='not run',
            blobtools_status='not run',whole_contig_candidates='',priority_intervals='',evidence_intervals='',
            source_bp='',protein_searched_bp='',megablast_searched_bp='',sensitive_blastn_status='not run')
        links = []
        if name in fcs:
            item = fcs[name]
            current = fcs_rows(item['actions'])
            counts = Counter(x['action'] for x in current)
            review_count = counts['REVIEW']+counts['REVIEW_RARE']
            actionable = sum(counts[k] for k in ('EXCLUDE','TRIM','FIX','REVIEW','REVIEW_RARE'))
            row.update(fcs_status='completed',fcs_actions=len(current),fcs_review_actions=review_count,
                       audit_status='No actionable GX findings' if not actionable else
                       ('REVIEW / CLEANUP REQUIRED' if final else 'See cleaning and review actions'))
            actions.extend(dict(assembly=name,coordinate_frame='FCS input; 1-based inclusive',**x) for x in current)
            for key, title in [('actions','GX actions'),('taxonomy','GX taxonomy'),('contaminants','GX removed sequence')]:
                if item.get(key):
                    links.append(f'[{title}](fcs/{quote(Path(item[key]).name)})')
        if name in adapters:
            current = adaptor_rows(adapters[name]['report'])
            count = sum(x.get('action') in ('ACTION_TRIM','ACTION_EXCLUDE') for x in current)
            row.update(adaptor_status='completed',adaptor_actions=count)
            if final and count:
                row['audit_status']='REVIEW / CLEANUP REQUIRED'
            elif final and not count and row['audit_status']=='No actionable GX findings':
                row['audit_status']='No actionable FCS findings'
            adapter_actions.extend(dict(assembly=name,coordinate_frame='adaptor input; 1-based inclusive',**x) for x in current)
            links.append(f"[Adaptor actions](fcs/{quote(Path(adapters[name]['report']).name)})")
        if name in metrics:
            current = rows(metrics[name]['metrics'])
            if any(x.get('classification_source') != 'FCS-GX' for x in current):
                raise ValueError('Plot evidence must explicitly identify FCS-GX as its source')
            fcs_sequences.extend(dict(cohort_assembly=name,**x) for x in current)
            links.append(f"[GC/depth and FCS assignments](fcs/{quote(Path(metrics[name]['metrics']).name)})")
        summary.append(row)
        lines.append(f"| {name} | {row['adaptor_actions'] if name in adapters else 'not run'} | "
                     f"{row['fcs_actions'] if name in fcs else 'not run'} | "
                     f"{row['fcs_review_actions'] if name in fcs else 'not run'} | {row['audit_status']} | "
                     f"{' · '.join(links) or 'Not run'} |")
        if name in blob:
            review = Path(blob[name]['review'])
            data = json.loads((review/'evidence.json').read_text())
            if data['assembly'] != name:
                raise ValueError('BlobTools evidence identity mismatch: '+name)
            current = rows(review/'review.tsv')
            intervals.extend(dict(cohort_assembly=name,**x) for x in current)
            current_sequences=rows(review/'sequence_evidence.tsv')
            sequences.extend(dict(cohort_assembly=name,**x) for x in current_sequences)
            row.update(blobtools_status='completed',whole_contig_candidates=data['whole_contig_candidates'],
                priority_intervals=sum(bool(x.get('support')) and x.get('priority') in
                    ('concordant non-host','extended non-host') for x in current),
                evidence_intervals=sum(bool(x.get('support')) for x in current),
                source_bp=data.get('source_bp',sum(int(x['length_bp']) for x in current_sequences)),
                protein_searched_bp=data['protein_searched_bp'],megablast_searched_bp=data['megablast_searched_bp'],
                sensitive_blastn_status=data['nucleotide_status'])
    lines += ['', '## GC and read coverage across assemblies', '']
    if fcs_sequences:
        lines += ['![GC and read coverage coloured by FCS-GX](report/cohort_gc_coverage.png)', '',
            'One panel per assembly and one dot per contig; dot area represents contig length. Colours use **FCS-GX '
            'divisions and explicit unresolved/mixed categories**, not BLAST or BlobTools assignments. Repeat, '
            'low-coverage and other inconclusive FCS results are not promoted to host assignments. '
            'All contigs are included where GC and read depth are measurable; zero-depth contigs are included. '
            'Panel labels count observations omitted for missing measurements. '
            '[PDF](report/cohort_gc_coverage.pdf) · [Plot scope](report/cohort_gc_coverage_scope.tsv) · '
            '[All FCS plot evidence](report/cohort_fcs_sequence_evidence.tsv)', '']
    else:
        lines += ['FCS GC/coverage plot unavailable: FCS plot evidence was not supplied.', '']
    if blob:
        lines += ['## Optional BlobTools evidence', '',
            'These are separately enabled, exhaustive similarity-search results. Higher-priority intervals have '
            'concordant or extended non-host evidence under the existing review rules. All signal intervals also '
            'include short signals and competing hits. Counts are not automatic removal decisions. '
            'Search columns give megabases and percent of the post-FCS assembly examined, not job counts or hit counts. '
            'The methods can cover the same bases, so their amounts must not be added.', '',
            '| Assembly | Whole-contig candidates | Higher-priority intervals | All signal intervals | Protein BLASTX searched: Mb (%) | MEGABLAST searched: Mb (%) | Results |',
            '|---|---:|---:|---:|---:|---:|---|']
        for row in summary:
            def amount(key):
                bp=int(row[key]); total=int(row['source_bp'])
                return f'{bp/1e6:.2f} ({100*bp/total:.2f}%)' if total else 'NA'
            name=row['assembly']; prefix=f'blobtools/{quote(name)}/'
            lines += [f"| {name} | {row['whole_contig_candidates']} | {row['priority_intervals']} | "
                f"{row['evidence_intervals']} | {amount('protein_searched_bp')} | {amount('megablast_searched_bp')} | "
                f"[Evidence report]({prefix}review/README.md) · [Decision sheet]({prefix}review/review.tsv) |"]
        lines += ['', 'Use the linked decision sheets to make explicit selections and preserve their source checksum, '
                  'assembly, coordinate stage and review ID. Sensitive BLASTN status is recorded in the cohort summary. '
                  'Historical sampled results, if inspected locally, remain sampled evidence and must not be described '
                  'as exhaustive. Production heuristic sampling has been retired.', '']
    lines += ['## Tables and follow-up', '',
        '[Assembly summaries](report/cohort_summary.tsv) · [GX actions](report/fcs_actions.tsv) · '
        '[Adaptor actions](report/adaptor_actions.tsv)', '',
        'Review GX REVIEW/REVIEW_RARE entries and any final-audit cleaning recommendations. Keep adaptor-input, '
        'GX-input and post-FCS coordinates separate. Apply approved changes in the matching upstream coordinate frame, '
        'then resume assembly and re-screen the resulting final FASTA. Do not silently trim finalized sequences while '
        'retaining old indexes, name maps, alignments or coordinate-based reports.', '']
    if blob:
        lines += ['[All BlobTools decision rows](report/cohort_review.tsv) · '
                  '[All BlobTools sequence evidence](report/cohort_sequence_evidence.tsv)', '']
    write_table(detail/'cohort_summary.tsv',summary,list(summary[0]) if summary else ['assembly','stage'])
    for filename, records, fallback in [
        ('fcs_actions.tsv',actions,['assembly','coordinate_frame','seq_id','start_pos','end_pos','action']),
        ('adaptor_actions.tsv',adapter_actions,['assembly','coordinate_frame','accession','action','range']),
        ('cohort_fcs_sequence_evidence.tsv',fcs_sequences,['cohort_assembly','sequence']),
        ('cohort_review.tsv',intervals,['cohort_assembly','action_id']),
        ('cohort_sequence_evidence.tsv',sequences,['cohort_assembly','sequence'])]:
        write_table(detail/filename,records,list(dict.fromkeys(k for x in records for k in x)) or fallback)
    (detail/'report.json').write_text(json.dumps(dict(stage=spec['stage'],fcs_enabled=spec['fcs_enabled'],
        blobtools_enabled=spec['blobtools_enabled'],plot_classification_source='FCS-GX' if fcs_sequences else 'unavailable',
        assemblies=summary),indent=2)+'\n')
    (out/'contamination_report.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    build(json.loads(args.spec.read_text()),args.out)
