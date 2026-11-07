"""Prepare full production PSMC comparisons from a successful preset-mask audit."""
import argparse
import csv
import json
from pathlib import Path
from prepare_integration_params import prepare


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();source=a.source.resolve();out=a.out.resolve();repo=Path(__file__).resolve().parents[2]
    m=json.loads((source/'manifest.json').read_text())
    status=json.loads((source/'reports/status.json').read_text())
    if m.get('mode')!='alignment_presets' or status['status']!='PASS':raise ValueError('Expected successful preset comparison')
    if out.exists() and any(out.iterdir()):raise ValueError('Use a new full-comparison directory')
    baseline=Path(m['baseline']);eligibility=baseline/'assembly_eligibility.json'
    data=json.loads(eligibility.read_text());sample=m['sample'];j=m['controls'][0]
    rows=sorted((r for r in data['assemblies'] if r['sample']==sample),key=lambda r:r['haplotype'])
    if len(rows)!=2 or [r['fasta'] for r in rows]!=j['assemblies']:raise ValueError('Baseline final assemblies differ from audited inputs')
    for r in rows:
        for key in ['fasta','fai']:
            if not Path(r[key]).is_file():raise ValueError('Missing final input: '+r[key])
    out.mkdir(parents=True,exist_ok=True)
    manifest=out/'individual.tsv'
    with manifest.open('w',newline='') as f:
        fields=['sample','mode','include_bed','exclude_bed']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        w.writerow(dict(sample=sample,mode='assemblies',include_bed=j.get('include_bed') or '',exclude_bed=j.get('exclude_bed') or ''))
    for preset in ['asm10','asm20']:
        dest=out/preset;dest.mkdir()
        params=prepare(baseline,repo/'nextflow.config',dest/'params.json')
        params.update({'psmc_'+k:v for k,v in j['filters'].items()})
        params.update(run_psmc=True,psmc_mode='assemblies',psmc_assembly_preset=preset,
            psmc_manifest=str(manifest),psmc_final_eligibility=str(eligibility),psmc_bootstraps=10,
            psmc_mutation_rate=None,psmc_generation_time=None,outdir=str(dest/'results'),publish_dir_mode='symlink')
        (dest/'params.json').write_text(json.dumps(params,indent=2)+'\n')
    (out/'comparison.json').write_text(json.dumps(dict(repository=str(repo),source=str(source),
        sample=sample,baseline=str(baseline),presets=['asm10','asm20']),indent=2)+'\n')
    print('Prepared full assembly-route comparisons for '+sample)


if __name__=='__main__':main()
