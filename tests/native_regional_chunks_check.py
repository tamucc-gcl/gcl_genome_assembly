"""Known-answer whole-CLI test of supported chunk settings, with no biological names."""
from pathlib import Path
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from pangenome_regional_patch import cli_regions, coordinate_audit

original,fixed,root=map(Path,sys.argv[1:]);root.mkdir(parents=True,exist_ok=True)
paf=root/'input.paf';query='synthetic_source';length=8000000
blocks=[(i*50000,i*50000+40000,'A') for i in range(60)]+[(3000000,length,'B')]
paf.write_text(''.join(f'{query}\t{length}\t{lo}\t{hi}\t+\t{chrom}\t10000000\t0\t{hi-lo}\t{hi-lo}\t{hi-lo}\t60\n' for lo,hi,chrom in blocks))
for variant,binary in [('native',original),('fixed',fixed)]:
    for chunk in (1000000,100000,25000):
        dest=root/f'{variant}_{chunk}';dest.mkdir()
        cmd=[str(binary.resolve()),'-p',str(paf.resolve()),'-b',str((dest/'split_').resolve()),
             '-n','.25','-Q','2','-u',str(chunk),'-s','-a','_AMBIGUOUS_','-L',str((dest/'log').resolve())]
        with (dest/'stderr').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=log)
        rows=cli_regions(dest,{query:length});coordinate_audit(rows,[dict(assembly='fixture',sequence=query,source_bp=length)])
        secondary=[r for r in rows if r['chromosome']=='A']
        if chunk!=25000:assert not secondary,(variant,chunk,secondary)
        elif variant=='fixed':
            assert len(secondary)==1,secondary
            assert secondary[0]['start_0based']==0 and secondary[0]['end_exclusive']==2990000,secondary
        else:assert len(secondary)==60,len(secondary)
print('PASS: real CLI chunk thresholds distinguish long runs of short intervals; grouping and ownership verified in six cases')
