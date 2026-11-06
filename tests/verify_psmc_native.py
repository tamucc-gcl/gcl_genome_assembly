"""Verify published native fixtures and scaling against upstream export tables."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import os
import shutil

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from psmc_results import export,read_fasta
from run_provenance import collect


def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--iteration',type=int,default=2)
    a=p.parse_args();published=a.root/'results/post_assembly/psmc'
    rows=list(csv.DictReader((published/'summary.tsv').open(),delimiter='\t'))
    assert {r['sample']:r['status'] for r in rows}=={'assembly_ind':'completed','read_ind':'completed','haploid_ind':'skipped'}
    for name in ['combined_history.png','combined_history.pdf','individual_histories.png',
                 'individual_histories.pdf','combined_curves.tsv']:
        assert (published/name).stat().st_size>0,name
    hap1=dict(read_fasta(a.root/'fixtures/hap1.fa'));hap2=dict(read_fasta(a.root/'fixtures/hap2.fa'))
    true_hets=sum(x!=y for x,y in zip(hap1['long_scaffold'],hap2['long_scaffold']))
    references=[]
    for sample in ['assembly_ind','read_ind']:
        d=published/'1'/sample
        for name in ['input/callable.bed','input/diploid.fa.gz','input/input.psmcfa','input/input_qc.json',
                     'input/provenance.json','report/history.png','report/history.pdf','report/report.md',
                     'report/curves.tsv','report/interval_qc.tsv','fits/observed.psmc']:
            assert (d/name).stat().st_size>0,name
        qc=json.loads((d/'input/input_qc.json').read_text())
        reference=json.loads((d/'reference/reference.json').read_text())
        expected=a.root/'fixtures'/('hap1.fa' if sample=='assembly_ind' else 'read_reference.fa')
        assert reference['sample']==sample and reference['id'].startswith(sample+'_')
        assert reference['sha256']==hashlib.sha256(expected.read_bytes()).hexdigest()==qc['reference_sha256']
        references.append(reference['sha256'])
        assert qc['status']=='ready' and qc['heterozygous_snps']>20 and qc['psmc_retained_callable_bp']>10000
        if sample=='assembly_ind':assert qc['heterozygous_snps']==true_hets
        else:assert .95*true_hets<=qc['heterozygous_snps']<=1.05*true_hets
        assert 'short_piece' not in (d/'input/callable.bed').read_text()
        if sample=='read_ind':assert qc['genotype_confidence_counts']['PL_separation']>0
        assert (d/'input/mask_audit.tsv').stat().st_size>0
        assert (d/'input/conversion_audit.tsv').stat().st_size>0
        seeds=list((d/'fits').glob('bootstrap_*.resampling.json'))
        fits=list((d/'fits').glob('bootstrap_*.psmc'))
        assert len(seeds)==len(fits)>0
        assert len({json.loads(p.read_text())['derived_seed'] for p in seeds})==len(seeds)
        native=a.root/(sample+'-upstream')
        command=['psmc_plot.pl','-R','-n'+str(a.iteration),'-s50','-u1e-8','-g0.5',str(native),str(d/'fits/observed.psmc')]
        # A gnuplot executable is unnecessary for numerical exports; it is not
        # claimed as a successful upstream figure render when unavailable.
        result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (a.root/(sample+'-upstream.log')).write_text(result.stdout)
        expected=[list(map(float,line.split())) for line in Path(str(native)+'.0.txt').read_text().splitlines() if line and not line.startswith('#')]
        actual=export([d/'fits/observed.psmc'],a.root/(sample+'-adapter-comparison'),a.iteration,50,1e-8,.5)
        assert len(expected)==len(actual),(len(expected),len(actual))
        for upstream,ours in zip(expected,actual):
            assert abs(upstream[0]-ours['time'])<=max(1e-5,abs(ours['time'])*1e-5)
            assert abs(upstream[1]*10000-ours['population'])<=max(1e-3,abs(ours['population'])*1e-5)
    assert len(set(references))==2,'Distinct individual references were collapsed'
    originals=sorted((a.root/'work').rglob('versions.tsv'))
    assert originals
    before=Path.cwd()
    with tempfile.TemporaryDirectory() as tmp:
        temp=Path(tmp);specs=[]
        for i,source in enumerate(originals,1):
            dest=temp/'staged'/('version'+str(i));dest.mkdir(parents=True)
            shutil.copyfile(source,dest/'versions.tsv');specs.append(dict(path=str(source)))
        try:
            os.chdir(temp);collect(specs,temp/'staged',dict(catalog_sources=[]),temp/'collected')
            issues=list(csv.DictReader((temp/'collected/software_version_issues.tsv').open(),delimiter='\t'))
            assert not issues,issues
        finally:os.chdir(before)
    print('PASS: both native modes, scope/QC, published figures, upstream scaling/iteration/terminal-interval equivalence and strict native version collection')


if __name__=='__main__':main()
