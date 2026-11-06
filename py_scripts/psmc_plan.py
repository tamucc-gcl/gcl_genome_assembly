"""Plan optional PSMC jobs using final biological identities and shared references."""
import argparse
import csv
import json
import math
from pathlib import Path
import re


def positive(value, name):
    if value in (None, '', 'null', 'None'):
        return None
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(name + ' must be positive and finite')
    return result


def manifest_rows(path):
    if not path or Path(path).name.startswith('NO_'):
        return {}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        allowed = {'sample', 'mode', 'include_bed', 'exclude_bed', 'mutation_rate', 'generation_time'}
        if not reader.fieldnames or 'sample' not in reader.fieldnames or set(reader.fieldnames)-allowed:
            raise ValueError('PSMC manifest requires sample and only documented columns')
        rows = {}
        for row in reader:
            sample = row['sample'].strip()
            if not sample or sample in rows:
                raise ValueError('Empty or duplicate PSMC manifest sample')
            row = {k: (v or '').strip() for k, v in row.items()}
            for key in ['include_bed', 'exclude_bed']:
                if row.get(key):
                    p = Path(row[key])
                    if not p.is_absolute():
                        p = Path(path).resolve().parent / p
                    row[key] = str(p.resolve(strict=True))
            rows[sample] = row
    return rows


def plan(eligibility, references, reads, overrides=None, mode='auto', mutation_rate=None, generation_time=None,min_contig_bp=1000000):
    overrides = overrides or {}
    refs = {str(r['taxid']): r for r in references}
    if len(refs) != len(references):
        raise ValueError('Duplicate shared-reference taxid')
    readmap = {}
    for r in reads:
        if r['sample'] in readmap:
            raise ValueError('Duplicate complete HiFi read input')
        readmap[r['sample']] = r['fastq']
    by_sample = {}
    for row in eligibility['assemblies']:
        by_sample.setdefault(row['sample'], []).append(row)
    if set(overrides)-set(by_sample):
        raise ValueError('Unknown PSMC manifest sample(s): '+', '.join(sorted(set(overrides)-set(by_sample))))
    selected = sorted(overrides if overrides else by_sample)
    if (mutation_rate is not None or generation_time is not None) and len({str(r['taxid']) for s in selected for r in by_sample[s]}) > 1:
        raise ValueError('For multiple species, supply scaling per sample in the PSMC manifest')
    jobs, skipped = [], []
    for sample in selected:
        rows = by_sample[sample]
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', sample):
            raise ValueError('Unsafe PSMC sample identifier: '+sample)
        taxa = {str(r['taxid']) for r in rows}
        if len(taxa) != 1 or len({r['ploidy'] for r in rows}) != 1:
            raise ValueError('Inconsistent within-individual PSMC identity: '+sample)
        taxid = next(iter(taxa))
        o = overrides.get(sample, {})
        chosen = o.get('mode') or mode
        if chosen not in ['auto', 'assemblies', 'hifi']:
            raise ValueError('PSMC mode must be auto, assemblies or hifi')
        mu = positive(o.get('mutation_rate') or mutation_rate, 'mutation_rate')
        gen = positive(o.get('generation_time') or generation_time, 'generation_time')
        if (mu is None) != (gen is None):
            raise ValueError('Supply both mutation rate and generation time, or leave both null')
        reason = None
        final = sorted((r for r in rows if r['completion'] == 'finalized'), key=lambda r:r['haplotype'])
        pair = len(final) == 2 and {r['haplotype'] for r in final} == {'hap1', 'hap2'}
        if rows[0]['ploidy'] != 2:
            reason = 'Standard PSMC requires a biologically diploid individual'
        elif len(final) != len(rows):
            reason = 'Expected finalized assembly representation is incomplete'
        elif taxid not in refs:
            reason = 'No shared final linear reference available'
        elif 'sequences' in refs[taxid] and not any(r['length']>=min_contig_bp for r in refs[taxid]['sequences']):
            reason = 'No reference scaffolds meet the PSMC minimum length'
        elif chosen == 'assemblies' and not pair:
            reason = 'Assembly mode requires both same-individual haplotypes'
        elif chosen == 'hifi' and sample not in readmap:
            reason = 'Read mode requires the complete same-individual HiFi read set'
        elif chosen == 'auto':
            if pair:
                chosen = 'assemblies'
            elif len(final) == 1 and sample in readmap:
                chosen = 'hifi'
            else:
                reason = 'Need two haplotypes, or one assembly plus complete HiFi reads; shotgun-only PSMC is unsupported'
        if reason:
            skipped.append(dict(sample=sample, taxid=taxid, status='skipped', reason=reason))
            continue
        ref = refs[taxid]
        if ref.get('coordinate_stage') != 'final' or not re.fullmatch(r'[0-9a-f]{64}', ref.get('sha256','')):
            raise ValueError('Shared PSMC reference lacks verified final identity')
        jobs.append(dict(id=sample, sample=sample, taxid=taxid, species=rows[0].get('species',''),
                         mode=chosen, reference=ref, assemblies=[r['fasta'] for r in final],
                         reads=readmap.get(sample), representation=final[0]['representation'],
                         include_bed=o.get('include_bed') or None, exclude_bed=o.get('exclude_bed') or None,
                         mutation_rate=mu, generation_time=gen,
                         scope_note='Explicit reference-coordinate inclusion BED' if o.get('include_bed') else
                         'Long reference scaffolds; autosomal status not established'))
    return dict(schema_version=1, jobs=jobs, skipped=skipped)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['eligibility','references','reads','out']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--manifest');p.add_argument('--mode',default='auto')
    p.add_argument('--mutation-rate');p.add_argument('--generation-time')
    p.add_argument('--min-contig-bp',type=int,default=1000000)
    a=p.parse_args()
    result=plan(json.loads(Path(a.eligibility).read_text()),json.loads(Path(a.references).read_text()),
                json.loads(Path(a.reads).read_text()),manifest_rows(a.manifest),a.mode,
                positive(a.mutation_rate,'mutation_rate'),positive(a.generation_time,'generation_time'),a.min_contig_bp)
    Path(a.out).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    Path('psmc_plan.md').write_text('# PSMC input plan\n\n'+str(len(result['jobs']))+' eligible individuals; '+
        str(len(result['skipped']))+' skipped.\n\n[Input and reference identities](psmc_plan.json).\n',encoding='utf-8')


if __name__=='__main__':
    main()
