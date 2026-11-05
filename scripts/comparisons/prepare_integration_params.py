"""Reuse recorded, currently supported settings for a cached integration check."""
import argparse
import json
from pathlib import Path
import re


def prepare(baseline,config,out):
    source=Path(baseline)/'provenance/run_provenance.json'
    data=json.loads(source.read_text(encoding='utf-8'))['context']['parameters']
    # Public defaults live in the first params block. Exclude process/profile settings.
    text=Path(config).read_text(encoding='utf-8').split('params {',1)[1].split('\n}',1)[0]
    supported=set(re.findall(r'^    (\w+)\s*=',text,re.M))
    # The caller supplies these invocation-specific inputs and destinations.
    supplied={'sample_sheet','hic_readsets','outdir','chimera_break','chimera_hifi_bam_manifest'}
    values={k:v for k,v in data.items() if k in supported-supplied}
    Path(out).write_text(json.dumps(values,indent=2)+'\n',encoding='utf-8')
    return values


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline','config','out']:p.add_argument('--'+name,required=True)
    a=p.parse_args();prepare(a.baseline,a.config,a.out)
