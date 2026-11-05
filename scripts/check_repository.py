#!/usr/bin/env python3
"""Check maintained source/document links and the frozen default-parameter inventory."""
import argparse
import csv
import io
from pathlib import Path
import re
from urllib.parse import unquote,urlsplit

ROOT=Path(__file__).resolve().parents[1]


def defaults():
    text=(ROOT/'nextflow.config').read_text(encoding='utf-8').split('params {',1)[1].split('\n}',1)[0]
    return [(m[1],m[2].split(' //',1)[0].strip()) for m in
            re.finditer(r'^    (\w+)\s*=\s*(.+)$',text,re.M)]


def inventory_text():
    buffer=io.StringIO();writer=csv.writer(buffer,delimiter='\t',lineterminator='\n')
    writer.writerow(['parameter','default_expression']);writer.writerows(defaults())
    return buffer.getvalue()


def check():
    errors=[];includes=0;links=0
    for folder in ['modules','workflows','functions','.']:
        for path in ([ROOT/'main.nf'] if folder=='.' else (ROOT/folder).glob('*.nf')):
            text=path.read_text(encoding='utf-8')
            text=re.sub(r'/\*.*?\*/|(?m:^\s*//[^\n]*$)','',text,flags=re.S)
            for relative in re.findall(r'\bfrom\s+[\'\"]([^\'\"]+)[\'\"]',text):
                target=path.parent/(relative if relative.endswith('.nf') else relative+'.nf')
                includes+=1
                if not target.is_file():errors.append(f'{path.relative_to(ROOT)}: missing include {relative}')
    for path in [ROOT/'main.nf',*(ROOT/'workflows').glob('*.nf'),*(ROOT/'modules').glob('*.nf')]:
        text=path.read_text(encoding='utf-8')
        text=re.sub(r'/\*.*?\*/|(?m:^\s*//[^\n]*$)','',text,flags=re.S)
        for relative in re.findall(r'\$\{projectDir\}/([^\'\"\n]+)',text):
            if '$' in relative or any(c in relative for c in '*?{}'):continue
            if Path(relative).suffix not in {'.py','.R','.nf','.yaml','.yml'}:continue
            if not (ROOT/relative).is_file():errors.append(f'{path.relative_to(ROOT)}: missing staged source {relative}')
    for path in [ROOT/'README.md',*(ROOT/'docs').rglob('*.md'),*(ROOT/'planning_docs').glob('*.md'),ROOT/'scripts/comparisons/README.md']:
        text=re.sub(r'```.*?```','',path.read_text(encoding='utf-8'),flags=re.S)
        for value in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',text):
            target=urlsplit(value.strip('<>'))
            if target.scheme or target.netloc or not target.path:continue
            links+=1
            if not (path.parent/unquote(target.path)).exists():errors.append(f'{path.relative_to(ROOT)}: missing link {value}')
    inventory=ROOT/'docs/parameter-defaults.tsv'
    if not inventory.is_file() or inventory.read_text(encoding='utf-8')!=inventory_text():
        errors.append('Parameter defaults changed; review the interface and explicitly refresh docs/parameter-defaults.tsv')
    return errors,includes,links


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-defaults',action='store_true');a=p.parse_args()
    if a.write_defaults:(ROOT/'docs/parameter-defaults.tsv').write_text(inventory_text(),encoding='utf-8',newline='\n')
    errors,includes,links=check()
    for error in errors:print(error)
    print(f'{"FAIL" if errors else "PASS"}: {includes} includes, {links} documentation links, {len(defaults())} configured defaults')
    raise SystemExit(bool(errors))
