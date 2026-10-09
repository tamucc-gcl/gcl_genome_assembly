#!/usr/bin/env python3
"""Explicit tool doubles: test interfaces and arithmetic, never mapping accuracy."""
import sys,re,time,subprocess
from pathlib import Path
tool=Path(sys.argv[0]).name;args=sys.argv[1:]
if args==['--version']:
    print(tool+' integration-double');sys.exit()


def fasta(path):
    records={};name=None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):name=line[1:].split()[0];records[name]=''
        else:records[name]+=line.strip()
    return records


if tool=='date':
    if args==['+%s%3N']:print(time.time_ns()//1000000)
    else:sys.exit(subprocess.call(['/usr/bin/date',*args]))
elif tool=='samtools' and args[0]=='faidx':
    path=Path(args[1]);records=fasta(path)
    if len(args)==2:Path(str(path)+'.fai').write_text(''.join(f'{n}\t{len(s)}\t0\t0\t0\n' for n,s in records.items()))
    else:
        n,region=args[2].rsplit(':',1);lo,hi=map(int,region.split('-'));print('>'+args[2]);print(records[n][lo-1:hi])
elif tool=='samtools' and args[0]=='sort':Path(args[args.index('-o')+1]).write_text(sys.stdin.read())
elif tool=='samtools' and args[0]=='index':Path(args[1]+'.bai').touch()
elif tool=='samtools' and args[0]=='quickcheck':assert Path(args[1]).read_text().startswith('@SQ')
elif tool=='samtools' and args[0]=='idxstats':
    for line in Path(args[1]).read_text().splitlines():
        if line.startswith('@SQ'):
            f=line.split('\t');print(f[1][3:]+'\t'+f[2][3:]+'\t2\t0')
elif tool=='samtools' and args[0]=='view':print(Path(args[2]).read_text(),end='')
elif tool=='minimap2' and '-ax' in args:
    for name,seq in fasta(args[-2]).items():
        print(f'@SQ\tSN:{name}\tLN:{len(seq)}')
        cigar=''.join(str(len(s))+('D' if s[0]=='N' else 'M') for s in re.findall('N+|[^N]+',seq))
        for i in range(2):print(f'r{i}\t0\t{name}\t1\t60\t{cigar}\t*\t0\t0\t*\t*\tNM:i:{seq.count("N")}')
elif tool=='minimap2':
    for q,seq in fasta(args[-1]).items():
        fragment=max(re.findall('[ACGT]+',seq),key=len)
        for t,target in fasta(args[-2]).items():
            start=target.find(fragment)
            if start>=0:
                lo=seq.find(fragment);n=len(fragment)
                print(f'{q}\t{len(seq)}\t{lo}\t{lo+n}\t+\t{t}\t{len(target)}\t{start}\t{start+n}\t{n}\t{n}\t60\tcg:Z:{n}=')
else:raise RuntimeError('Unexpected tool call: '+str(sys.argv))
