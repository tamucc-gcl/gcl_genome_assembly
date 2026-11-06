"""Label observed versions independently of PSMC analysis and presentation helpers."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['versions'])
    p.add_argument('--tools',nargs='+',required=True)
    p.add_argument('--out',default='versions.tsv')
    a=p.parse_args()
    rows=[]
    for tool in a.tools:
        if tool=='dipcall':
            records=list((Path(sys.prefix)/'conda-meta').glob('dipcall-*.json'))
            if len(records)!=1:raise ValueError('Cannot identify installed Dipcall package')
            version=json.loads(records[0].read_text())['version']+' (installed conda package)'
        elif tool=='python':version=sys.version.split()[0]
        else:
            cmd={'psmc':['psmc'],'minimap2':['minimap2','--version'],
                 'samtools':['samtools','--version-only'],'bcftools':['bcftools','--version-only']}[tool]
            run=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if tool=='psmc':
                matches=re.findall(r'^Version:\s*(\S+)',run.stdout,re.M)
                if len(matches)!=1:raise ValueError('Cannot identify PSMC version')
                version=matches[0]
            else:
                if run.returncode or not run.stdout.strip():raise ValueError('Version command failed: '+tool)
                version=run.stdout.splitlines()[0]
        rows.append(tool+'\t'+version)
    Path(a.out).write_text('\n'.join(rows)+'\n')


if __name__=='__main__':main()
